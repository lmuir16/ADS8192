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

# Create a GitHub repository and push
use_github()

# This will:
# 1. Create a new repo on GitHub
# 2. Add the remote
# 3. Push your initial commit

# If this doesn't work, you can create the repo manually on GitHub
# and then add the remote via the terminal:
# git remote add origin https://github.com/your-username/ADS8192.git
# git push -u origin main

# MIT License is a good default for open source
use_mit_license()

# This adds:
# - LICENSE.md file
# - License field in DESCRIPTION
