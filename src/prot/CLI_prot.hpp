#pragma once
#include "config_prot.hpp"
#include <string>

// Parse argc/argv into ProtConfig.
// Throws std::invalid_argument on bad input.
// Prints usage to stderr and returns false if -h/--help was given.
bool parse_prot_cli(int argc, char** argv, ProtConfig& cfg);
void print_prot_usage(const char* prog);
