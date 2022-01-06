if (interactive()) {
  suppressMessages(require(devtools))
  suppressMessages(require(usethis))
  suppressMessages(require(testthat))
  suppressMessages(require(remotes))
  
  suppressMessages(credentials::set_github_pat())
  
  # add git branch to console prompt
  if (requireNamespace("prompt", quietly = TRUE)) {
    prompt_git <- function(...) {
      paste0(
        "[", prompt::git_branch(), "]",
        " > "
      )
    }
    prompt::set_prompt(prompt_git)
    rm(prompt_git)
  }
}

options(
  tidyverse.quiet = TRUE
)
