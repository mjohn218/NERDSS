Example simulation input files for flat clathrin sheet formation in solution

- `parms_rev.inp` features reversible ring closure (default), where association rates of species within the same complex are calculated as a function of copy numbers in the complex
- `parms_irrev.inp` features irreversible ring closure (denoted by reaction block keyword `irrevRingClosure`), which sets association rates of species within the same complex to unity
- Several rates used:
  - K<sub>D</sub> = 100:
    - k<sub>a</sub> = 0.0332 nm<sup>3</sup>/μs
    - k<sub>b</sub> = 1.00022 s<sup>-1</sup>
  - K<sub>D</sub> = 1:
    - k<sub>a</sub> = 3.3971 nm<sup>3</sup>/μs
    - k<sub>b</sub> = 1.0225 s<sup>-1</sup>
  - K<sub>D</sub> = 0.2:
    - k<sub>a</sub> = 18.67 nm<sup>3</sup>/μs
    - k<sub>b</sub> = 1.124 s<sup>-1</sup>
