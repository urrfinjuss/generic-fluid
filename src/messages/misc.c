#include "gfluid.h"
#include "config.h"

void display_help() {
  printf("run ./gfluid with exactly one (1) command line argument: filename\n");
  printf("where the file contains configuration of the simulation.\n");
}

void display_version() {
  printf("%s\n", PACKAGE_STRING);
}

