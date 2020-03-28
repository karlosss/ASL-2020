#include <bits/stdc++.h>
#include "modes.h"

int main(int argc, char** argv) {
    if(argc >= 2 && strcmp(argv[1], "test") == 0) test_mode();
    else if(argc >= 2 && strcmp(argv[1], "debug") == 0) debug_mode();
    else if(argc >= 2 && strcmp(argv[1], "runtime") == 0) runtime_mode();
    else if(argc >= 2 && strcmp(argv[1], "flops") == 0) flops_mode();
}
