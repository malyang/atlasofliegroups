# auto_cluster_algebra research deliveries

This directory is a persistent fallback delivery channel for the Hall-algebra / quantum-cluster-algebra research project.

The requested standalone repository `malyang/auto_cluster_algebra` cannot be created by the currently installed GitHub App because the connector exposes repository-content writes but not the account-level create-repository permission.  Until that repository is created once in the GitHub UI, this dedicated branch publishes stable files here.

## Stable download

- `releases/latest.zip` — newest TeX + PDF + code research package
- `releases/SHA256SUMS.txt` — checksums
- `source/` — source used by the automatic build

Raw stable URL after the workflow finishes:

```text
https://raw.githubusercontent.com/malyang/atlasofliegroups/auto_cluster_algebra-artifacts/auto_cluster_algebra/releases/latest.zip
```

The GitHub Actions workflow on this branch compiles the paper, executes the exact symbolic verification, assembles a versioned ZIP, and commits both the versioned ZIP and `latest.zip` back to this branch.
