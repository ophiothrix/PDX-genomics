library(ProjectTemplate)
create.project("../PDX-genomics/", merge.strategy = "allow.non.conflict")

## Clean auxilary R files
system("rm ./diagnostics/1.R ./lib/helpers.R ./munge/01-A.R ./profiling/1.R ./src/eda.R ./tests/1.R")
system("rm TODO ./cache/README.md ./config/README.md ./data/README.md ./diagnostics/README.md ./doc/README.md ./graphs/README.md ./lib/README.md ./logs/README.md ./munge/README.md ./profiling/README.md ./reports/README.md ./src/README.md ./tests/README.md")
