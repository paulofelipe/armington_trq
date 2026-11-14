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

Whether it actually reproduces those results has not yet been verified.

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
uv run python src/armington_trq.py
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

- `src/armington_trq.py` – the `PathMCP` solver and the `ArmingtonTRQ` model.

## Known assumptions

- Region 0 is treated as the domestic origin in the welfare calculation, and
  its initial tariff `t0[0]` is assumed to be 0.
- Quota regions must satisfy `t_out >= t_in`.
