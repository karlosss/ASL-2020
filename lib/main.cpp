#include <bits/stdc++.h>
#include "modes.h"

int main(int argc, char** argv) {
    if(argc >= 2 && strcmp(argv[1], "test") == 0) test_mode();
    else if(argc >= 2 && strcmp(argv[1], "execute") == 0) run_mode();
    else if(argc >= 2 && strcmp(argv[1], "flops") == 0) flops_mode();
}
