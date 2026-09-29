# Donor random-intercept mixed models

These scripts retain the uploaded models, data handling, condition labels, contrasts, and multiplicity settings. Paths are configurable; the HbF spike file header was corrected and its existing result tables are now saved as CSV.

All three fit `log2(response) ~ condition + (1 | donor)` and calculate fold change as `2^estimate`. This is the ratio corresponding to a contrast on the log2 scale, not a ratio of raw arithmetic means. Each script fits its own dataset; do not combine their results into a single model.

| Script | Default input / response | Control | Comparisons and adjustment |
| --- | --- | --- | --- |
| run_um_model.R | day3_postediting.csv / value_day3 | unedited mock | Selected target conditions only; explicit adjust="dunnett" |
| run_mixed_model.R | day14_long_format.csv / value_day14 | CCR5_CBE | All conditions versus control using trt.vs.ctrl's default adjustment; four specified tEPOR contrasts with Holm |
| run_mixed_model_spike.R | HbF_spike.csv / HbF | Ctrl | All conditions versus Ctrl with explicit adjust="dunnettx"; four specified tEPOR contrasts with Holm |

The exact adjustment arguments above are retained, including their differences. The generic Dunnett wording in the source comments should not be used to claim identical computation across scripts. The day14 script relies on the installed emmeans default for trt.vs.ctrl; record the package version and resulting output before describing the precise implementation in Methods. The source spike script loads lme4 and emmeans; the other two additionally load lmerTest. Package loading and degrees-of-freedom settings were not standardized because that could change results.

## Inputs and execution

CSV files must be in long format with `donor`, `condition`, and the relevant numeric response column. Condition names must exactly match the settings in each script. Edit the response/column settings for other datasets only intentionally.

```r
Sys.setenv(MODEL_INPUT = "/path/to/day14_long_format.csv",
           MODEL_OUTPUT = "/path/to/results/day14")
source("run_mixed_model.R")
```

Or from Terminal, with the relevant R packages installed:

```bash
MODEL_INPUT=/path/to/HbF_spike.csv MODEL_OUTPUT=/path/to/results/spike Rscript run_mixed_model_spike.R
```

Run each dataset in a fresh R session for consistent package loading. Requires lme4 and emmeans, plus lmerTest for run_um_model.R and run_mixed_model.R. Record the actual package/R versions used in the study.

## Original data handling retained

- run_um_model.R removes missing donor/condition/response rows, removes nonfinite or nonpositive responses, skips absent targets, and fits only the control plus available requested targets. Its post-contrast label filter is retained; check the output contains the requested comparisons.
- run_mixed_model.R removes nonfinite or nonpositive responses; original factor handling and missing donor/condition behavior are retained.
- run_mixed_model_spike.R stops on missing or nonpositive HbF values; it retains the original fitting behavior for missing donor/condition and other invalid values.

No new data imputation, replicate averaging, tests, or model terms were added. Holm adjustment is performed across the four specified within-group contrasts separately from the control-comparison family.

## Outputs

run_um_model.R writes mixedmodel_vs_unedited_mock.csv.

run_mixed_model.R writes mixedmodel_vs_control.csv and mixedmodel_paired_contrasts.csv.

run_mixed_model_spike.R now writes HbF_spike_vs_Ctrl.csv and HbF_spike_within_group_contrasts.csv using the same tables it already computed and printed. This output-only addition does not modify the analyses.

## Verification limits

The statistical bodies were compared with the original scripts after accounting for path, header, and CSV-output edits. R and the analysis CSVs were unavailable in the preparation environment, so parsing, model execution, convergence, and numerical agreement could not be tested here.
