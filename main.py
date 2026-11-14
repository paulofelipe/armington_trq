"""Run the policy scenarios of Hallren & Riker (2018), "A Second Comparison of
Partial Equilibrium Models of TRQs with Sensitivity Analysis" (USITC Economics
Working Paper 2018-02-B), for the endogenous-price Armington model.

Scenarios 1-4 (Tables 1, 3-6) are run and compared with the paper's
"Endogenous Price Armington" column. Scenario 5 (Tables 7-8) is not run: its
reported ranges (e.g. a 0% domestic price change) are not those of the
endogenous-price Armington model.

Usage:  uv run python main.py [--plot]
"""

import sys

sys.path.insert(0, "src")

from armington_trq import ArmingtonTRQ

# Table 2 inputs are the ArmingtonTRQ defaults: shares 60/30/10, market size
# 100, eps = (1, 10, 10), sigma = 4, eta = -1. Only the policy changes by scenario.
# Regions: domestic, subject (faces the TRQ), non-subject.
SCENARIOS = {
    1: dict(
        name="Standard 10% ad valorem tariff",
        qbar=(0, 0, 0),
        t_in=(0, 0.10, 0),
        t_out=(0, 0.10, 0),
    ),
    2: dict(
        name="TRQ binding at 10% (in 0%, out 10%, quota 15)",
        qbar=(0, 15, 0),
        t_in=(0, 0.00, 0),
        t_out=(0, 0.10, 0),
    ),
    3: dict(
        name="TRQ binding at the quota (in 10%, out 50%, quota 15)",
        qbar=(0, 15, 0),
        t_in=(0, 0.10, 0),
        t_out=(0, 0.50, 0),
    ),
    4: dict(
        name="Non-binding TRQ (in 40%, out 50%, quota 15)",
        qbar=(0, 15, 0),
        t_in=(0, 0.40, 0),
        t_out=(0, 0.50, 0),
    ),
}

# Endogenous Price Armington column of Tables 3-6, as published (rounded).
PAPER = {
    1: dict(
        price_dom_pct=2.1,
        price_subj_pct=7.9,
        price_overall_pct=3.6,
        volume_subj=24.6,
        output_dom_pct=2.1,
        volume_subj_pct=-17.9,
    ),
    2: dict(
        price_dom_pct=2.1,
        price_subj_pct=7.9,
        price_overall_pct=3.6,
        volume_subj=24.6,
        output_dom_pct=2.1,
        volume_subj_pct=-17.9,
    ),
    3: dict(
        price_dom_pct=6.4,
        price_subj_pct=28.5,
        price_overall_pct=10.9,
        volume_subj=15.0,
        output_dom_pct=6.4,
        volume_subj_pct=-50.0,
    ),
    4: dict(
        price_dom_pct=6.7,
        price_subj_pct=30.1,
        price_overall_pct=11.4,
        volume_subj=14.5,
        output_dom_pct=6.7,
        volume_subj_pct=-51.8,
    ),
}

LABELS = {
    "price_dom_pct": "Pct chg in price of dom production (%)",
    "price_subj_pct": "Pct chg in price of subj imports (%)",
    "price_overall_pct": "Pct chg in overall prices (%)",
    "volume_subj": "Volume of subject imports",
    "output_dom_pct": "Pct chg in domestic production (%)",
    "volume_subj_pct": "Pct chg in volume of subj imports (%)",
}


def run_scenario(k):
    cfg = SCENARIOS[k]
    model = ArmingtonTRQ(**{key: v for key, v in cfg.items() if key != "name"})
    return model.solve()


def print_comparison(k, model):
    results, paper = model.summary(), PAPER[k]
    print(f"\nScenario {k}: {SCENARIOS[k]['name']}")
    print(f"{'':<42}{'This code':>10}{'Paper':>9}{'Diff':>8}")
    for key, label in LABELS.items():
        print(
            f"{label:<42}{results[key]:>10.1f}{paper[key]:>9.1f}"
            f"{results[key] - paper[key]:>8.1f}"
        )
    print(f"PATH status: {model.info['status']}")


if __name__ == "__main__":
    for k in SCENARIOS:
        model = run_scenario(k)
        print_comparison(k, model)
        if "--plot" in sys.argv:
            model.plot_welfare(save=f"welfare_scenario_{k}.png")
