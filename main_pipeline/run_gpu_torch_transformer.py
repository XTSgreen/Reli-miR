import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from sklearn.decomposition import PCA


class PCTransformer(torch.nn.Module):
    def __init__(
        self,
        seq_len: int,
        output_dim: int,
        d_model: int,
        nhead: int,
        num_layers: int,
        dim_feedforward: int,
        dropout: float,
    ) -> None:
        super().__init__()
        self.token_proj = torch.nn.Linear(1, d_model)
        self.pos_embed = torch.nn.Parameter(torch.zeros(1, seq_len, d_model))
        encoder_layer = torch.nn.TransformerEncoderLayer(
            d_model=d_model,
            nhead=nhead,
            dim_feedforward=dim_feedforward,
            dropout=dropout,
            activation="gelu",
            batch_first=True,
            norm_first=True,
        )
        self.encoder = torch.nn.TransformerEncoder(encoder_layer, num_layers=num_layers)
        self.head = torch.nn.Sequential(
            torch.nn.LayerNorm(d_model),
            torch.nn.Linear(d_model, d_model),
            torch.nn.GELU(),
            torch.nn.Dropout(dropout),
            torch.nn.Linear(d_model, output_dim),
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        x = x.unsqueeze(-1)
        x = self.token_proj(x) + self.pos_embed
        x = self.encoder(x)
        x = x.mean(dim=1)
        return self.head(x)


def load_matrix(path: Path) -> tuple[np.ndarray, np.ndarray, list[str]]:
    df = pd.read_csv(path)
    sample_ids = df.iloc[:, 0].astype(str).to_numpy()
    feature_names = df.columns[1:].tolist()
    values = df.iloc[:, 1:].to_numpy(dtype=np.float32)
    return sample_ids, values, feature_names


def mean_r2(y_true: np.ndarray, y_pred: np.ndarray) -> np.ndarray:
    rss = np.sum((y_true - y_pred) ** 2, axis=0)
    centered = y_true - y_true.mean(axis=0, keepdims=True)
    tss = np.sum(centered ** 2, axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        r2 = 1.0 - rss / tss
    r2[~np.isfinite(r2)] = np.nan
    return r2


def pearson_by_target(y_true: np.ndarray, y_pred: np.ndarray) -> np.ndarray:
    out = np.full(y_true.shape[1], np.nan, dtype=np.float32)
    for j in range(y_true.shape[1]):
        a = y_true[:, j]
        b = y_pred[:, j]
        if np.std(a) > 1e-8 and np.std(b) > 1e-8:
            out[j] = np.corrcoef(a, b)[0, 1]
    return out


def rmse_by_target(y_true: np.ndarray, y_pred: np.ndarray) -> np.ndarray:
    return np.sqrt(np.mean((y_true - y_pred) ** 2, axis=0))


def build_metric_frame(
    dataset: str,
    model: str,
    target_names: list[str],
    y_true: np.ndarray,
    y_pred: np.ndarray,
) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "dataset": dataset,
            "model": model,
            "miRNA": target_names,
            "R2": mean_r2(y_true, y_pred),
            "Pearson": pearson_by_target(y_true, y_pred),
            "RMSE": rmse_by_target(y_true, y_pred),
        }
    )


def build_global_frame(metric_df: pd.DataFrame) -> pd.DataFrame:
    group = metric_df.groupby("model", as_index=False).agg(
        mean_R2=("R2", "mean"),
        median_R2=("R2", "median"),
        mean_Pearson=("Pearson", "mean"),
        median_Pearson=("Pearson", "median"),
        mean_RMSE=("RMSE", "mean"),
        median_RMSE=("RMSE", "median"),
        n_miRNA=("miRNA", "count"),
    )
    return group.sort_values(["median_Pearson", "mean_R2", "mean_RMSE"], ascending=[False, False, True])


def to_device(x: np.ndarray, device: torch.device) -> torch.Tensor:
    return torch.tensor(x, dtype=torch.float32, device=device)


def save_prediction_csv(path: Path, sample_ids: np.ndarray, target_names: list[str], preds: np.ndarray) -> None:
    pd.DataFrame({"patient_id": sample_ids, **{name: preds[:, i] for i, name in enumerate(target_names)}}).to_csv(
        path,
        index=False,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-dir", default=".", help="StrongMiRRecon 项目目录")
    parser.add_argument("--model-name", default="Torch_GPU_Transformer", help="模型名称")
    parser.add_argument("--epochs", type=int, default=400)
    parser.add_argument("--patience", type=int, default=40)
    parser.add_argument("--pca-dim", type=int, default=128)
    parser.add_argument("--d-model", type=int, default=64)
    parser.add_argument("--nhead", type=int, default=8)
    parser.add_argument("--num-layers", type=int, default=2)
    parser.add_argument("--dim-feedforward", type=int, default=256)
    parser.add_argument("--dropout", type=float, default=0.2)
    parser.add_argument("--lr", type=float, default=1e-3)
    parser.add_argument("--weight-decay", type=float, default=1e-4)
    args = parser.parse_args()

    project_dir = Path(args.project_dir).resolve()
    input_dir = project_dir / "output" / "gpu_model_inputs"
    output_dir = project_dir / "output" / "gpu_transformer"
    output_dir.mkdir(parents=True, exist_ok=True)

    train_ids, x_train, feature_names = load_matrix(input_dir / "train_x_z.csv")
    val_ids, x_val, _ = load_matrix(input_dir / "validation_x_z.csv")
    test_ids, x_test, _ = load_matrix(input_dir / "test_x_z.csv")

    _, y_train, target_names = load_matrix(input_dir / "train_y_z.csv")
    _, y_val, _ = load_matrix(input_dir / "validation_y_z.csv")
    _, y_test, _ = load_matrix(input_dir / "test_y_z.csv")

    n_components = min(args.pca_dim, x_train.shape[0] - 1, x_train.shape[1])
    pca = PCA(n_components=n_components, svd_solver="full", random_state=123)
    x_train_pca = pca.fit_transform(x_train).astype(np.float32)
    x_val_pca = pca.transform(x_val).astype(np.float32)
    x_test_pca = pca.transform(x_test).astype(np.float32)

    if args.d_model % args.nhead != 0:
        raise ValueError("d-model 必须能被 nhead 整除。")

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    model = PCTransformer(
        seq_len=x_train_pca.shape[1],
        output_dim=y_train.shape[1],
        d_model=args.d_model,
        nhead=args.nhead,
        num_layers=args.num_layers,
        dim_feedforward=args.dim_feedforward,
        dropout=args.dropout,
    ).to(device)

    optimizer = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)
    loss_fn = torch.nn.MSELoss()

    train_x_tensor = to_device(x_train_pca, device)
    train_y_tensor = to_device(y_train, device)
    val_x_tensor = to_device(x_val_pca, device)
    val_y_tensor = to_device(y_val, device)
    test_x_tensor = to_device(x_test_pca, device)

    best_val = float("inf")
    best_epoch = 0
    best_state = None
    history: list[dict[str, float]] = []

    for epoch in range(1, args.epochs + 1):
        model.train()
        optimizer.zero_grad(set_to_none=True)
        pred_train = model(train_x_tensor)
        train_loss = loss_fn(pred_train, train_y_tensor)
        train_loss.backward()
        optimizer.step()

        model.eval()
        with torch.no_grad():
            pred_val = model(val_x_tensor)
            val_loss = loss_fn(pred_val, val_y_tensor).item()

        history.append(
            {
                "epoch": epoch,
                "train_loss": float(train_loss.item()),
                "validation_loss": float(val_loss),
            }
        )

        if val_loss < best_val:
            best_val = val_loss
            best_epoch = epoch
            best_state = {k: v.detach().cpu().clone() for k, v in model.state_dict().items()}
        elif epoch - best_epoch >= args.patience:
            break

    if best_state is None:
        raise RuntimeError("Transformer 模型没有得到有效状态。")

    model.load_state_dict(best_state)
    model.eval()

    with torch.no_grad():
        pred_val = model(val_x_tensor).detach().cpu().numpy()
        pred_test = model(test_x_tensor).detach().cpu().numpy()

    metric_val = build_metric_frame("validation", args.model_name, target_names, y_val, pred_val)
    metric_test = build_metric_frame("test", args.model_name, target_names, y_test, pred_test)

    global_val = build_global_frame(metric_val)
    global_test = build_global_frame(metric_test)

    save_prediction_csv(output_dir / f"{args.model_name}_prediction_validation_z.csv", val_ids, target_names, pred_val)
    save_prediction_csv(output_dir / f"{args.model_name}_prediction_test_z.csv", test_ids, target_names, pred_test)

    metric_val.to_csv(output_dir / f"{args.model_name}_validation_metrics_by_miRNA.csv", index=False)
    metric_test.to_csv(output_dir / f"{args.model_name}_test_metrics_by_miRNA.csv", index=False)
    global_val.to_csv(output_dir / f"{args.model_name}_validation_global_metrics.csv", index=False)
    global_test.to_csv(output_dir / f"{args.model_name}_test_global_metrics.csv", index=False)
    pd.DataFrame(history).to_csv(output_dir / f"{args.model_name}_training_history.csv", index=False)

    torch.save(
        {
            "model_state_dict": best_state,
            "pca_mean": pca.mean_,
            "pca_components": pca.components_,
            "pca_explained_variance_ratio": pca.explained_variance_ratio_,
            "feature_names": feature_names,
            "target_names": target_names,
            "config": vars(args),
            "best_epoch": best_epoch,
            "best_validation_loss": best_val,
        },
        output_dir / f"{args.model_name}.pt",
    )

    base_val = pd.read_csv(project_dir / "output" / "11_validation_global_metrics.csv")
    base_test = pd.read_csv(project_dir / "output" / "13_test_global_metrics.csv")

    combined_val = pd.concat([base_val, global_val], ignore_index=True).sort_values(
        ["median_Pearson", "mean_R2", "mean_RMSE"],
        ascending=[False, False, True],
    )
    combined_test = pd.concat([base_test, global_test], ignore_index=True).sort_values(
        ["median_Pearson", "mean_R2", "mean_RMSE"],
        ascending=[False, False, True],
    )

    combined_val.to_csv(output_dir / "combined_validation_global_metrics_with_transformer.csv", index=False)
    combined_test.to_csv(output_dir / "combined_test_global_metrics_with_transformer.csv", index=False)

    summary = {
        "model_name": args.model_name,
        "device": str(device),
        "best_epoch": best_epoch,
        "best_validation_loss": best_val,
        "pca_dim": int(n_components),
        "d_model": int(args.d_model),
        "nhead": int(args.nhead),
        "num_layers": int(args.num_layers),
        "dim_feedforward": int(args.dim_feedforward),
        "dropout": float(args.dropout),
        "lr": float(args.lr),
        "weight_decay": float(args.weight_decay),
    }
    (output_dir / f"{args.model_name}_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print(json.dumps(summary, ensure_ascii=False))


if __name__ == "__main__":
    main()
