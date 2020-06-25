# Configuration

Before running a new project, the default configuration should be adjusted. To do this, create a new config file named "config.yml" in the root folder of the project. Alternatively, config files can be provided by running snakemake with the `--config <file>` command.

The project configuration file overwrites the default configs found in `config/config_defaults.yml`. If you want to run only assembly for example, it should be enough to create a new config file containing:

```
steps:
    - assembly
```

