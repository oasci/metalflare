import json
import os
import string

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from deeptime.markov.msm import MaximumLikelihoodMSM

from metalflare.analysis.figures import use_mpl_rc_params

os.chdir(os.path.dirname(os.path.realpath(__file__)))

# Specify the paths to the trajectory and topology files
base_dir = "../../../"

# Update plot params
rc_json_path = os.path.join(base_dir, "misc/003-figure-style/matplotlib-rc-params.json")
font_dirs = [os.path.join(base_dir, "misc/003-figure-style/roboto")]
use_mpl_rc_params(rc_json_path, font_dirs)  # type: ignore

# === Input files === #
label_paths = {
    "Reduced": os.path.join(
        base_dir, "analysis/009-pw-configs/data/reduced-clustering-labels.npy"
    ),
    "Oxidized": os.path.join(
        base_dir, "analysis/009-pw-configs/data/oxidized-clustering-labels.npy"
    ),
    "Cu(I)": os.path.join(
        base_dir, "analysis/009-pw-configs/data/cu-clustering-labels.npy"
    ),
}

output_dir = "./"
os.makedirs(output_dir, exist_ok=True)

# === Settings === #
lagtime = 49
min_transition_prob = 0.001
min_occ_prob = 0.001

min_edge_width = 0.7
edge_width_scale = 5.0
min_head_width = 8
head_width_scale = 12

# Each simulation states has three concatenated trajectories that are the same
# length, but with different initial conditions.
# We need to mark these trajectories boundaries to avoid any nonphysical
# transitions.
n_independent_runs = 3

# Load cluster centers with labels (A, B, C, ...)
cluster_json_path = os.path.join(
    base_dir,
    "analysis/009-pw-configs/data/vamp-clustering-global.json",
)
with open(cluster_json_path, "r") as f:
    cluster_centers = json.load(f)["cluster_centers"]


# === Create subplots === #
fig, axes = plt.subplots(1, 3, figsize=(12, 4))

colors = {"Reduced": "#1e2e79", "Oxidized": "#EC4067", "Cu(I)": "#f99752"}


def get_occ(state, occ, state_symbols):
    if state in state_symbols:
        i_occ = np.argwhere(state_symbols == state)[0][0]
        val = occ[i_occ]
    else:
        val = 0.0
    return val


def get_trans_prob(state_i, state_j, trans_mat, state_symbols):
    if state_i in state_symbols:
        i = np.argwhere(state_symbols == state_i)[0][0]
    else:
        return 0.0
    if state_j in state_symbols:
        j = np.argwhere(state_symbols == state_j)[0][0]
    else:
        return 0.0

    return trans_mat[i, j]


n_states = len(cluster_centers)
for ax, (name, path) in zip(axes, label_paths.items()):
    traj = np.load(path)
    state_letters = list(string.ascii_uppercase[:n_states])

    # Split by run
    n_steps_per_traj = int(traj.shape[0] / n_independent_runs)
    run_idxs = np.repeat(np.arange(n_independent_runs), n_steps_per_traj)
    traj_runs = np.split(traj, n_independent_runs)

    msm = MaximumLikelihoodMSM(reversible=True).fit_fetch(traj_runs, lagtime=lagtime)
    trans_prob = msm.transition_matrix
    occ = msm.stationary_distribution

    G = nx.DiGraph()
    state_symbols = msm.state_symbols()
    for i in range(n_states):
        size = get_occ(i, occ, state_symbols)
        G.add_node(i, size=size, label=state_letters[i])

    for i in range(n_states):
        occ_i = get_occ(i, occ, state_symbols)
        if occ_i == 0.0:
            continue
        for j in range(n_states):
            occ_j = get_occ(j, occ, state_symbols)
            if occ_j == 0.0:
                continue
            if i == j:
                continue
            p = get_trans_prob(i, j, trans_prob, state_symbols)
            if (
                (p > min_transition_prob)
                and (occ_i > min_occ_prob)
                and (occ_j > min_occ_prob)
            ):
                G.add_edge(i, j, weight=p)

    pos = nx.circular_layout(G)
    min_node_size = 200
    node_sizes = {n: max(4000 * G.nodes[n]["size"], min_node_size) for n in G.nodes}

    nx.draw_networkx_nodes(
        G,
        pos,
        ax=ax,
        node_size=node_sizes.values(),
        node_color=colors[name],
        edgecolors="k",
    )
    labels = {n: G.nodes[n]["label"] for n in G.nodes}
    nx.draw_networkx_labels(
        G,
        pos,
        labels=labels,
        font_size=10,
        font_weight="bold",
        ax=ax,
        font_color="#ffffff",
    )

    for u, v in G.edges():
        weight = G[u][v]["weight"]
        scaled_weight = np.sqrt(weight)  # Linear scaling did not look correct visually
        width = max(min_edge_width, edge_width_scale * scaled_weight)
        head_width = min_head_width + head_width_scale * scaled_weight

        target_radius = (node_sizes[v] / np.pi) ** 0.5
        min_target_margin = target_radius * 0.85

        rad = 0.3 if G.has_edge(v, u) and u != v else -0.2

        nx.draw_networkx_edges(
            G,
            pos,
            edgelist=[(u, v)],
            width=[width],
            arrows=True,
            arrowstyle="->",
            connectionstyle=f"arc3,rad={rad}",
            min_source_margin=0,
            min_target_margin=min_target_margin,
            arrowsize=head_width,
            ax=ax,
        )

    ax.set_title(name, color="black")
    ax.axis("off")
    ax.set_xlim(-1.3, 1.2)
    ax.set_ylim(-1.3, 1.2)

plt.tight_layout()
plt.subplots_adjust(left=0.00, right=0.99, top=0.99, bottom=0.05)

y_label = 0.53
x_delta = 0.02
fig.text(
    x=0.14 + x_delta,
    y=y_label,
    s="Reduced",
    ha="center",
    fontsize=12,
    fontweight="bold",
)
fig.text(
    x=0.475 + x_delta,
    y=y_label,
    s="Oxidized",
    ha="center",
    fontsize=12,
    fontweight="bold",
)
fig.text(
    x=0.81 + x_delta, y=y_label, s="Cu(I)", ha="center", fontsize=12, fontweight="bold"
)

out_path = os.path.join(output_dir, "fig006.svg")
plt.savefig(out_path)
plt.close()
