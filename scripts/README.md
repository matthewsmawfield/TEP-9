# TEP-EFA Analysis Scripts

This folder contains the reproducible analysis pipeline for Paper 15 (Earth Flyby).

## Structure

```
scripts/
├── steps/           # Numbered analysis steps
│   ├── step_01_*.py # Data acquisition
│   ├── step_02_*.py # Data processing
│   └── ...
└── utils/           # Shared utilities
```

## Execution

Steps should be run in numerical order. Each step produces outputs in `results/outputs/`.
