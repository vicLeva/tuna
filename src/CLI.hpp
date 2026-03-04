#pragma once

#include "Config.hpp"

void print_usage(const char* prog);

// Parses argv into cfg.  Returns false (and prints an error) on bad input;
// returning false is also used for --help so the caller should show usage.
bool parse_args(int argc, char* argv[], Config& cfg);
