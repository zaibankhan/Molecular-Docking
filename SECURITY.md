# Security & Responsible Use

Thank you for taking the time to help keep this project safe and used responsibly.

## 1. Reporting a software vulnerability

If you discover a security issue in the code (e.g. remote code execution in the
web server, credential/data leakage, unsafe file handling), **do not** open a
public issue. Instead, report it privately.

- Prefer a [GitHub private security advisory](https://github.com/zaibankhan/Molecular-Docking/security/advisories/new)
  on the repository.
- Alternatively, open an issue with the `[SECURITY]` prefix and request that a
  maintainer triage it privately.

Please include:
- a concise description of the bug and its impact,
- the affected component/version,
- minimal steps to reproduce,
- any suggested fix (optional).

We aim to acknowledge reports within 5 working days.

## 2. Project scope and what this tool is NOT

This project is a **scientific/educational molecular docking pipeline**. It
prepares structures and runs the widely published, open-source **AutoDock Vina**
scoring function to rank candidate poses. It does **not**:

- make any toxicity, drug-safety, or bioactivity claims,
- infer that a predicted pose is biologically relevant,
- validate any compound for clinical, industrial, or consumer use.

Predicted scores are estimates for relative ranking in research/teaching, and
must not be treated as experimental measurements.

## 3. Responsible use

Molecular docking is a standard, broadly accessible academic technique with
legitimate uses in drug discovery, structural biology, and education. We ask
users to apply this tool and its results ethically and legally:

- **Do not** use predicted scores as the sole basis for decisions that affect
  human health or the environment (drug dosing, toxicity clearance, product
  formulation).
- **Respect** applicable laws and export/dual-use regulations for any compound
  chemistry you work with.
- **Do not** use the tool to facilitate the creation of chemical, biological, or
  nuclear weapons, or any harmful substance.
- **Attribute** realistically: predictions are computational estimates from the
  AutoDock Vina scoring function, not measured binding constants. Cite Vina when
  using results in research.
- **Publish transparently**: report the exact inputs, box, exhaustiveness, and
  seed so results are reproducible.

We believe in open science and will not restrict access to this standard tool,
but we ask all users to act responsibility and to comply with their local and
national laws.

## 4. Built-in ethical safeguards in the code

- Every report was designed to label affinities as **predicted**, not measured,
  and to list interpretation caveats.
- Interactive results pages in the web UI carry the same "predicted, not
  measured" language.
- The demo complex (streptavidin–biotin) is a textbook teaching example.

If you believe a feature encourages misuse, please report it under Section 1.