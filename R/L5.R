# Install development tools if needed
install.packages(c("devtools", "usethis", "roxygen2", "testthat"))

library(usethis)

# Create the package (choose your own name!)
# This creates a new directory with the package structure
create_package("~/ADS8192")  # Or wherever you want it

# This will open a new RStudio session in the package directory

# Initialize git repository
use_git()

# This will:
# 1. Create a .git directory
# 2. Make an initial commit
# 3. Potentially restart RStudio
