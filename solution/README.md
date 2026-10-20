# Reference approach

The reference implementation solves the three technical parts of the task:

1. The temperature dependence of the equilibrium constant is evaluated with the van't Hoff relation, and the reverse rate constant is derived from the forward rate constant and that equilibrium constant.
2. The fixed-temperature mass balance includes the reactor residence time and a nonzero B feed; the energy balance includes feed sensible heat and heat exchanged with the utility at `T - T_utility`.
3. The scan and bracket search retain every steady state in the requested temperature interval, and local stability is obtained from the two-state concentration/temperature Jacobian.

The corresponding checks are separated by scoring weight in `tests/test_outputs.py`: point checks cover the first part, heat-balance cases cover the second, and the multiplicity/stability cases cover the third.
