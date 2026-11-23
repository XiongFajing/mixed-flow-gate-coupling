# Gate-boundary coupling in partitioned mixed-flow models

MATLAB programs and input data for the numerical experiments in the associated manuscript.

## Requirements

MATLAB R2024a. No parallel execution is required.

## Running the examples

Open an example folder as the MATLAB current folder and run its entry point.

| Folder                     | Example                             | Entry point       |
| -------------------------- | ----------------------------------- | ----------------- |
| `01_single_pipe`           | Single-pipe benchmark               | `RUN_SINGLE`      |
| `02_two_pipe_explicit`     | Explicit two-pipe simulations       | `RUN_EXPLICIT`    |
| `03_two_pipe_implicit`     | Implicit two-pipe simulations       | `RUN_IMPLICIT`    |
| `04_two_pipe_convergence`  | Coupling convergence experiment     | `RUN_CONVERGENCE` |
| `05_smoothing_sensitivity` | Regularization-parameter experiment | `RUN_SMOOTHING`   |
| `06_canal_tunnel/explicit` | Explicit canal-tunnel simulation    | `RUN_CASE`        |
| `06_canal_tunnel/implicit` | Implicit canal-tunnel simulation    | `RUN_CASE`        |
| `07_local_feedback`        | Local feedback experiment           | `RUN_FEEDBACK`    |
| `08_gate_opening_cycle`    | Gate-opening variation experiment   | `RUN_OPENING`     |

For example, enter `01_single_pipe` and run:

```matlab
RUN_SINGLE
```

Restart MATLAB before switching between example folders to avoid conflicts between functions and classes with the same name. Do not add the entire repository recursively to the MATLAB path.

Parameters are specified in each example's configuration file or main simulation script.

## Outputs

Outputs are saved in each example's `results` folder as CSV and MAT files. Column headers identify variables and units. Water heads are relative to the local invert; gate flows apply over the intervals specified in the output tables. Check the status fields for simulation completion. Rerunning an example may overwrite existing outputs.

## Citation and contact

Please cite the associated manuscript when using this repository. Citation information is provided in `CITATION.cff`.

Contact: Fajing Xiong, 2018302060243@whu.edu.cn.
