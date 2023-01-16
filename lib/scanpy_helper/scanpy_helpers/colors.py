import altair as alt
import pandas as pd
from natsort import natsorted


def set_scale_anndata(adata, column, palette=None):
    if palette is None:
        palette = column

    adata._sanitize()

    tmp_cols = getattr(COLORS, palette)
    adata.uns[f"{column}_colors"] = [
        tmp_cols[cat] for cat in adata.obs[column].cat.categories
    ]


def altair_scale(variable, *, data=None, data_col=None, **kwargs):
    """
    Discrete color scale for altair based on our color definitions.

    Parameters:
    -----------
    variable
        name of the color scale
    data
        Data frame used for the chart. If specified, will only show values that actually occur in the data.
    data_col
        If specified, check this column in `data` instead of `variable`

    Returns
    -------
    Altair color scale
    """
    tmp_colors = getattr(COLORS, variable)
    if data is not None:
        data_col = variable if data_col is None else data_col
        tmp_colors = {k: tmp_colors[k] for k in sorted(data[data_col].unique())}

    return alt.Scale(
        domain=list(tmp_colors.keys()),
        range=list(tmp_colors.values()),
        **kwargs,
    )


def altair_scale_mpl(scheme, **kwargs):
    """
    Use a continuous color scheme from mpl with altair
    """
    from matplotlib import cm
    from matplotlib.colors import to_hex

    return alt.Scale(
        range=[to_hex(x) for x in cm.get_cmap(scheme, 1000)(range(1000))], **kwargs
    )


def plot_palette(variable):
    """Display a palette"""
    tmp_cols = getattr(COLORS, variable)
    return (
        alt.Chart(
            pd.DataFrame.from_dict(tmp_cols, orient="index", columns=["color"])
            .reset_index()
            .rename(columns={"index": variable})
        )
        .mark_rect(height=40, width=30)
        .encode(
            x=alt.X(variable),
            color=alt.Color(variable, scale=altair_scale(variable), legend=None),
        )
    )


def plot_all_palettes():
    return alt.vconcat(
        *[plot_palette(v) for v in COLORS.__dict__.keys() if not v.startswith("_")]
    ).resolve_scale(color="independent")

class COLORS:

    group = {
        "Ctrl": "#e41a1c",
        "A8": "#377eb8",
        "A9": "#4daf4a",
        "A8A9": "#984ea3",
    }

    sample = {
        "AJ1": "#a6cee3",
        "AJ2": "#1f78b4",
        "AJ3": "#b2df8a",
        "AJ4": "#33a02c",
        "AJ5": "#fb9a99",
        "AJ6": "#e31a1c",
        "AJ7": "#fdbf6f",
        "AJ8": "#ff7f00",
        "AJ9": "#cab2d6",
        "AJ10": "#6a3d9a",
        "AJ11": "#ffff99",
        "AJ12": "#b15928",
        "AJ13": "#fdb462",
        "AJ14": "#b3de69",
        "AJ15": "#ffed6f",
        "AJ16": "#ccebc5",
    }

    batch = {
        "1": "#ef8a62",
        "2": "#67a9cf",
    }

    sex = {
        "1": "#ef8a62",
        "2": "#67a9cf",
    }

    phase = {
        "G1": "#8da0cb",
        "G2M": "#66c2a5",
        "S": "#fc8d62",
    }

    cell_type = dict(
        natsorted(
            {
                "B cell": "#b5bd61",
                "B cell dividing": "#f7b6d2",
                "Enterocyte prox": "#8c564b",
                "Enterocyte dist": "#17becf",
                "Enteroendocrine": "#ffbb78",
                "TA Progenitor": "#aa40fc",
                "Goblet": "#1f77b4",
                "Paneth": "#e377c2",
                "Stem Progenitor": "#555555",
                "T cells CD4": "#279e68",
                "T cells CD8": "#006d2c",
                "Tuft": "#98df8a",
                "MO DC": "#aec7e8",
            }.items()
        )
    )

    cell_type_rough = dict(
        natsorted(
            {
                "B cell": "#b5bd61",
                "Progenitor": "#f7b6d2",
                "Enterocyte": "#8c564b",
                "Goblet": "#17becf",
                "Enteroendocrine": "#ffbb78",
                "Paneth": "#aa40fc",
                "Tuft": "#1f77b4",
                "MO DC": "#e377c2",
                "T cell": "#279e68",
            }.items()
        )
    )