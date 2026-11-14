# Armington TRQ

A partial-equilibrium Armington model with a tariff-rate quota (TRQ), solved as a
Mixed Complementarity Problem (MCP) with a PATH-style solver.

> **⚠️ AI experiment – not validated.**
> This repository is an experiment in AI-assisted code generation. The code has
> **not** been independently validated or checked against reference results.
> Do not rely on its output for research, policy or any other decision without
> verifying it yourself.

## Purpose

The code is intended to replicate the results in the USITC Economics Working
Paper:

<https://www.usitc.gov/publications/332/working_papers/ecwp-2018-02-b-02-21-18-as-pdf.pdf>

`main.py` runs scenarios 1–4 of the paper (Tables 1–6) and prints them next to the
paper's "Endogenous Price Armington" column; the values match to the one decimal
reported. Scenario 5 (Tables 7–8, sensitivity analysis) is not implemented: its
reported ranges are not those of the endogenous-price Armington model. The
exogenous-price Armington and DSK columns are also out of scope.
This is still a self-check by the same AI-generated code, not an independent
validation.

## Comparison with the paper

Output of `main.py` against the "Endogenous Price Armington" column of Tables 3–6
(Hallren & Riker, 2018), shown as "this code / paper".

| Metric | S1: 10% tariff | S2: in 0% / out 10% | S3: in 10% / out 50% | S4: in 40% / out 50% |
|---|---|---|---|---|
| | This code / Paper | This code / Paper | This code / Paper | This code / Paper |
| Price of dom. production (%) | 2.1 / 2.1 | 2.1 / 2.1 | 6.4 / 6.4 | 6.7 / 6.7 |
| Price of subject imports (%) | 7.9 / 7.9 | 7.9 / 7.9 | 28.5 / 28.5 | 30.1 / 30.1 |
| Overall prices (%) | 3.6 / 3.6 | 3.6 / 3.6 | 10.9 / 10.9 | 11.4 / 11.4 |
| Volume of subject imports | 24.6 / 24.6 | 24.6 / 24.6 | 15.0 / 15.0 | 14.5 / 14.5 |
| Domestic production (%) | 2.1 / 2.1 | 2.1 / 2.1 | 6.4 / 6.4 | 6.7 / 6.7 |
| Volume of subject imports (%) | -17.9 / -17.9 | -17.9 / -17.9 | -50.0 / -50.0 | -51.8 / -51.8 |

Scenario 2 (binding TRQ, 0% in-quota) equals scenario 1 in the Armington model,
as the paper notes: only the effective marginal tariff matters.

## What it does

- Total import demand is CES (Armington), with elasticity of substitution `sigma`
  across origins and aggregate demand elasticity `eta`.
- Each origin's supply responds to its received price with elasticity `eps`.
- Regions with a quota (`qbar > 0`) pay `t_in` in-quota and `t_out` out-of-quota.
  The quota shadow price (`omega`) is endogenous, bounded in `[0, t_out - t_in]`.
- Reports equilibrium price and quantity changes by region and a welfare
  decomposition for the domestic origin (consumers, producers, tariff revenue,
  quota rent), with a waterfall plot.

## Usage

Requires Python ≥ 3.12. With [uv](https://docs.astral.sh/uv/):

```bash
uv sync
uv run python main.py          # scenarios 1-4 vs. the paper
uv run python main.py --plot   # also save welfare_scenario_<k>.png
```

Or from Python:

```python
from armington_trq import ArmingtonTRQ

model = ArmingtonTRQ(
    regions=("dom", "sub", "nsub"),
    v0=(60.0, 30.0, 10.0),      # base trade values (quantities, since p0 = 1)
    eta=1.0,                     # total demand elasticity
    sigma=4.0,                   # Armington elasticity
    eps=(1.0, 10.0, 10.0),       # supply elasticities
    qbar=(0.0, 15.0, 0.0),       # quota volume (0 = no quota)
    t_in=(0.0, 0.1, 0.0),        # in-quota tariff
    t_out=(0.0, 0.5, 0.0),       # out-of-quota tariff
)
model.solve()
model.report()                   # tables printed to the terminal
model.plot_welfare()             # or plot_welfare(save="welfare.png")
```

## Layout

- `main.py` – runs the paper's scenarios and compares them with its tables.
- `src/armington_trq.py` – the `PathMCP` solver and the `ArmingtonTRQ` model.

## Known assumptions

- Region 0 is treated as the domestic origin in the welfare calculation, and
  its initial tariff `t0[0]` is assumed to be 0.
- Quota regions must satisfy `t_out >= t_in`.
