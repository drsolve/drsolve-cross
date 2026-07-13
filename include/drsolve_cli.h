#ifndef DRSOLVE_CLI_H
#define DRSOLVE_CLI_H

void drsolve_cli_print_version(void);
void drsolve_cli_print_usage(const char *prog_name);
int drsolve_cli_main(int argc, char *argv[], const char *prog_name);

#endif
