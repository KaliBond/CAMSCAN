import pytest

from cams_scoring_engine import (
    CAMS_NODES,
    NodeMetrics,
    check_divergence,
    compute_node_value,
    compute_pairwise_bond,
    compute_standard_bs,
    score_society,
)


def sample_nodes():
    return {
        "Helm": NodeMetrics(4, 5, 9, 5, bs_native=2),
        "Shield": NodeMetrics(3, 4, 10, 4, bs_native=1),
        "Lore": NodeMetrics(7, 7, 4, 8, bs_native=6),
        "Stewards": NodeMetrics(5, 6, 7, 5, bs_native=4),
        "Craft": NodeMetrics(6, 6, 5, 7, bs_native=5),
        "Hands": NodeMetrics(4, 5, 8, 4, bs_native=2),
        "Archive": NodeMetrics(7, 7, 4, 8, bs_native=6),
        "Flow": NodeMetrics(7, 7, 5, 7, bs_native=6),
    }


def test_node_value_formula():
    node = NodeMetrics(coherence=6, capacity=5, stress=4, abstraction=8)
    assert compute_node_value(node) == 11


def test_pairwise_bond_formula():
    a = NodeMetrics(6, 6, 4, 8)
    b = NodeMetrics(4, 5, 6, 5)
    # numerator=((6+4)*0.6)+((8+5)*0.4)=11.2 ; avg_stress=5 ; bond=11.2/6
    assert compute_pairwise_bond(a, b) == pytest.approx(11.2 / 6)


def test_standard_bs_formula():
    node = NodeMetrics(7, 6, 4, 5)
    assert compute_standard_bs(node) == 10


def test_check_divergence():
    assert check_divergence(2, 10) is True
    assert check_divergence(8, 10) is False
    assert check_divergence(None, 10) is False


def test_score_society_output_shape_and_columns():
    df = score_society("USA", 1861, sample_nodes())
    assert len(df) == 8
    assert df["Node"].tolist() == CAMS_NODES
    assert {
        "Society",
        "Year",
        "Node",
        "Coherence",
        "Capacity",
        "Stress",
        "Abstraction",
        "Node Value",
        "Bond Strength",
        "BS_Std",
        "BS_Native",
        "Divergence",
    }.issubset(df.columns)


def test_score_society_requires_all_nodes():
    nodes = sample_nodes()
    nodes.pop("Flow")
    with pytest.raises(ValueError, match="All 8 CAMS nodes must be provided"):
        score_society("USA", 1861, nodes)


def test_validate_enforces_bounds():
    bad_nodes = sample_nodes()
    bad_nodes["Helm"] = NodeMetrics(11, 5, 5, 5)
    with pytest.raises(ValueError, match="coherence must be between 1 and 10"):
        score_society("USA", 1861, bad_nodes)
