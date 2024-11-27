include(FetchContent)
FetchContent_Declare(
  lattice
  GIT_REPOSITORY https://github.com/todo-group/lattice.git
  GIT_TAG        master
)
list(APPEND FetchContents lattice)
