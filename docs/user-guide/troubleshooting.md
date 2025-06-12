# Troubleshooting

Recommended steps to troubleshoot the pipeline.

## 1.1 Email

Check your email for an email regarding pipeline failure. You will receive an email from slurm@biowulf.nih.gov with the subject: Slurm Job_id=[#] Name=CARLISLE Failed, Run time [time], FAILED, ExitCode 1

## 1.2 Review the log files

You can check logs in two ways to diagnose workflow errors:

Master SLURM Log File
Located in /path/to/results/dir/ and named slurm-[jobid].out, this file summarizes the overall workflow run. It identifies which rule failed and provides error details for any jobs executed locally via SLURM.

Individual Rule Log Files
After identifying the failed rule(s) from the master SLURM log, examine the corresponding log files in /path/to/results/dir/logs/.

Each file follows this naming convention:

{rulename}.{masterjobID}.{individualruleID}.{wildcards}.{out or err}`

.out files capture standard output
.err files capture standard error messages

## 1.3 Restart the run

Once you have identified and addressed the issue, you may resume the SINCLAIR run.

Unlock the output directory, perform another dry-run, and check the status of the pipeline, then resubmit to the cluster.

```
sinclair run \
    -profile biowulf \
    --input assets/input_manifest.csv \
    --contrast assets/contrast_manifest.csv \
    --output /data/$USER/scRNA_test \
    -params-file assets/params.yml
```

## 1.4 Help & Contributing

Come across a **bug**? Open an [issue](https://github.com/CCBR/SINCLAIR/issues) and include a minimal reproducible example.

Have a **question**? Ask it in [discussions](https://github.com/CCBR/SINCLAIR/discussions).

Want to **contribute** to this project? Check out the [contributing guidelines](../contributing.md).

**General Inquiries and Collaboration:** Please contact the CCBR Pipeliner team at [CCBR_Pipeliner@mail.nih.gov](mailto:CCBR_Pipeliner@mail.nih.gov).
