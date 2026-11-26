# Duality of Regularized System Identification

MATLAB code accompanying the paper **Duality of Time and Frequency Domain Kernel-Based Regularized System Identification**, by **Xiaozhu Fang, Yu Xu, and Tianshi Chen**. The testing databanks (1.6 GB) are saved in https://sandbox.zenodo.org/records/648679.  

The paper establishes sufficient conditions for duality between time-domain and frequency-domain kernel-based regularized estimates, including hyperparameter estimation. This repository provides discrete-time numerical examples using the diagonal/correlated (DC) kernel, with nonparametric transient modeling and comparisons with alternative transient treatments and a prediction error method (PEM).

Set `configs{3}` near the top of either experiment script to choose the number of MultiStart runs per hyperparameter. With five hyperparameters, the default `configs{3} = 5` gives 25 starts. Set it to `1` for five starts.

## Repository contents

| File or folder | Purpose |
| --- | --- |
| [`Automatica_main.m`](Automatica_main.m) | Original experiment driver: 2,500 systems in each of four data banks. Requires the full input dataset and local path/parallel-pool configuration. |
| [`Automatica_Dataresult/`](Auto_Data/) | Original saved results for the four data banks, including fit scores, objective values, and hyperparameters. These files can be inspected without rerunning the experiments. |
| [`example_for_small_data.m`](example_for_small_data.m) | Serial illustration using only 10 systems from S1D1*. |
| [`s1d1star_sysfrom1to10.mat`](s1d1star_sysfrom1to10.mat) | Input data and true models for the 10-system illustration. |
| [`rfir_eb_td.m`](rfir_eb_td.m) | Time-domain kernel-based regularized estimator. |
| [`rfir_eb_fd.m`](rfir_eb_fd.m) | Frequency-domain kernel-based regularized estimator. |
| [`rfir_eb_td_grad.m`](rfir_eb_td_grad.m), [`rfir_eb_fd_grad.m`](rfir_eb_fd_grad.m) | Time- and frequency-domain estimators with an analytic marginal-likelihood gradient and Hessian. |
| [`LICENSE`](LICENSE) | MIT license for the code. |

The `Dataset_TAC2025_` prefix in saved filenames is retained from the original implementation.

## Requirements

Use MATLAB with the following toolboxes:

- **System Identification Toolbox**: data/model objects, `n4sid`, `tfest`, and simulation.
- **Control System Toolbox**: impulse responses, frequency responses, and `lsim`.
- **Optimization Toolbox**: `fmincon`.
- **Global Optimization Toolbox**: `MultiStart`.
- **Parallel Computing Toolbox**: required by `Automatica_main.m` as written; the 10-system example runs serially.

## Quick start: 10-system illustration

Clone or download this repository:

```bash
git clone https://github.com/PI-DD-MCD-Lab/Duality-of-Time-and-Frequency-Domain-Kernel-Based-Regularized-System-Identification.git
```

In MATLAB, set the current folder to the repository root, then run:

```matlab
rng(0);  % Fix the random stream for this serial illustration.
example_for_small_data
```

The script loads the included `s1d1star_sysfrom1to10.mat`, processes 10 systems, prints iteration progress, and saves:

```text
Dataset_TAC2025_s1d1star_sysfrom1to10.mat
```

## Original experiments and full input data

The full experiments use four data banks, each containing 2,500 systems and corresponding input/output records. The sampling interval is 1, and the impulse-response fit is evaluated over the first 125 coefficients.

| Data bank | Required input file | System dynamics | Samples used per record | SNR |
| --- | --- | --- | --- | --- |
| S1D1* | `s1d1star.mat` | Fast | 500 | 10 |
| S1D2* | `s1d2star.mat` | Fast | 375 | 1 |
| S2D1* | `s2d1star.mat` | Slow | 500 | 10 |
| S2D2* | `s2d2star.mat` | Slow | 375 | 1 |

Each input file must contain `sys` and `dat` cell arrays: `sys{ii}` is the true system, and `dat{ii}` is its input/output record. The driver uses the last 500 or 375 samples of each record.

**Full input dataset download: pending.** A public download link or dataset DOI has not yet been added. The approximately 1.6 GB input dataset is required to rerun all four data banks; the files in `Auto_Data/` contain results, not the `sys` and `dat` inputs.

### Running the full experiment

1. Obtain the full input dataset and place the four input files in one folder. The following layout uses a local `raw_data/` folder:

   ```text
   repository-root/
   ├── Automatica_main.m
   ├── rfir_eb_td.m
   ├── rfir_eb_fd.m
   ├── Auto_Data/                  # Included original results
   └── raw_data/                   # Full input data, downloaded separately
       ├── s1d1star.mat
       ├── s1d2star.mat
       ├── s2d1star.mat
       └── s2d2star.mat
   ```

2. In `Automatica_main.m`, replace the author-specific `data_paths` assignment with the following. Keep the subsequent `valid_path` check and `cd` call. Adding the repository folder to the MATLAB path keeps the estimator functions accessible after the directory change.

   ```matlab
   repo_dir = fileparts(mfilename('fullpath'));
   addpath(repo_dir);
   data_paths = {fullfile(repo_dir, 'raw_data')};
   ```

   Alternatively, set `data_paths` to the absolute path of an external data folder.

3. Replace the original `parcluster` / `c.NumWorkers = 75` / `parpool(75)` block with a pool appropriate for your machine. For example, use the local profile's configured worker count:

   ```matlab
   if isempty(gcp('nocreate'))
       parpool('local');
   end
   ```

   To run serially, remove the pool setup and replace `parfor ii = 1:num_iters` with `for ii = 1:num_iters`.

4. Run `Automatica_main` from the repository root. It processes all four banks and saves `Dataset_TAC2025_<bank>.mat` in the selected data folder. A rerun writes to the same result filenames.


## Estimation configurations and outputs

The current scripts evaluate the following eight configurations in order. EIR denotes the extended impulse response model; FRF denotes the frequency response function model.

| Column | Domain | Configuration |
| --- | --- | --- |
| 1 | Frequency | FRF model with transient estimation |
| 2 | Frequency | EIR model with transient estimation |
| 3 | Time | FRF model with transient estimation |
| 4 | Time | EIR model with transient estimation |
| 5 | Time | EIR model without transient estimation (`eir-no-t`; prewindowing comparison) |
| 6 | Time | EIR model with truncation setting 125  （Not included in the paper）|
| 7 | Time | Truncated FIR model with the first 125 output samples discarded (`cut`) |
| 8 | Time | PEM baseline using `n4sid` followed by `tfest` |

Saved result files contain:

| Variable | Contents |
| --- | --- |
| `gFIT` | Impulse-response fit over the first 125 coefficients. |
| `yFIT` | Simulated-output fit over the last 125 samples of a newly generated 6,250-sample input sequence. |
| `GINF` | Maximum sampled FRF error normalized by an approximation of the true FRF's L2 norm. |
| `OBJ` | Hyperparameter optimization objective values. |
| `HP` | Hyperparameter vectors returned by the optimizers. |
| `algorithm` | Estimator configurations and their ordering. |



## Citation and license

If you use this code, please cite:

> Xiaozhu Fang, Yu Xu, and Tianshi Chen. *Duality of Time and Frequency Domain Kernel-Based Regularized System Identification*. Submitted to Automatica. 2026

The code is distributed under the [MIT License](LICENSE).
