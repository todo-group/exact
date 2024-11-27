include(FetchContent)
FetchContent_Declare(
  standards
  GIT_REPOSITORY https://github.com/todo-group/standards.git
  GIT_TAG        master
)
list(APPEND FetchContents standards)
