#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/sysinfo.h>
#include <sys/resource.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#define MAX_CPUS 64
unsigned long long prev_total_cpu_time[MAX_CPUS] = {0};
unsigned long long prev_idle_time[MAX_CPUS] = {0};

void get_cpu_usage(double *cpu_usage) {
    FILE *fp = fopen("/proc/stat", "r");
    if (fp == NULL) {
        perror("Failed to open /proc/stat");
        exit(EXIT_FAILURE);
    }

    char buffer[1024];
    char *line;
    int icpu=0;
    while ((line = fgets(buffer, sizeof(buffer), fp)) != NULL) {
        if (strncmp(line, "cpu", 3) == 0) {
            unsigned long long user, nice, system, idle, iowait, irq, softirq, steal, guest, guest_nice;
            sscanf(line, "cpu %llu %llu %llu %llu %llu %llu %llu %llu %llu %llu",
                   &user, &nice, &system, &idle, &iowait, &irq, &softirq, &steal, &guest, &guest_nice);

            unsigned long long total_cpu_time = user + nice + system + idle + iowait + irq + softirq + steal + guest + guest_nice;
            unsigned long long idle_time = idle + iowait;

	    icpu++;
            if (total_cpu_time == prev_total_cpu_time[icpu]) {
	      cpu_usage[icpu] = 0.001;
	    } else {
	      unsigned long long total_diff = total_cpu_time - prev_total_cpu_time[icpu];
	      unsigned long long idle_diff = idle_time - prev_idle_time[icpu];
	      cpu_usage[icpu] = ((total_diff - idle_diff) / (double)total_diff) * 100.0+1E-3;
	    }
            prev_total_cpu_time[icpu]=total_cpu_time;
            prev_idle_time[icpu]=idle_time;
	    
        }
    }

    fclose(fp);
}

int memory_usage()
{
  struct sysinfo sys_info;
  if(sysinfo(&sys_info) != 0) {
    perror("sysinfo");
    exit(EXIT_FAILURE);
  }
  
  long long total_memory = sys_info.totalram * sys_info.mem_unit;
  long long free_memory = sys_info.freeram * sys_info.mem_unit;
  long long available_memory = sys_info.freeram * sys_info.mem_unit + sys_info.bufferram * sys_info.mem_unit + sys_info.sharedram * sys_info.mem_unit;
  
  printf(" %8.3lf", (double) (total_memory)/(1024*1024*1024));
  printf(" %8.3lf", (double)(free_memory)/(1024*1024*1024));
  printf(" %8.3lf\n", (double)(available_memory)/(1024*1024*1024));
  
  return 0;
}

int main() {
    double cpu_usage[MAX_CPUS];
    long i=0;
    while (1) {
        get_cpu_usage(cpu_usage);

        printf("%08ld",i);
        for (int i = 0; i < MAX_CPUS; ++i) {
            if (cpu_usage[i] != 0.0) {
                printf(" %6.2f",cpu_usage[i]);
            }
        }
        memory_usage();
        sleep(1); // Wait for 1 second before checking CPU usage again
	i++;
    }

    return 0;
}


