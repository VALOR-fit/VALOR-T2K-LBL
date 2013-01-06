#! /usr/bin/env python


f = open('stats','r')

unknown = success = fail = running = scheduled = cleared = abort = cancelled = 0

for line in f:
    if "======================= glite-wms-job-status " in line:
        if not "Success =====================" in line:
            unknown = unknown + 1
            continue
    if "Current Status:" in line:
        if "Done (Success)" in line:
            success = success + 1
        elif "Running" in line:
            running = running + 1
        elif "Done (Failed)" in line:
            fail = fail + 1
        elif "Aborted" in line:
            abort = abort + 1
        elif "Scheduled" in line:
            scheduled = scheduled + 1
        elif "Cleared" in line:
            cleared = cleared + 1
        elif "Cancelled" in line:
            cancelled = cancelled + 1
        else:
            print line,
            unknown += 1

print unknown, " Unknown jobs"
print fail, " Failed jobs"
print abort, " Aborted jobs"
print cancelled, " Cancelled jobs"
print scheduled, " Scheduled jobs"
print running, " Running jobs"
print success, " Success jobs"
print cleared, " Cleared jobs"
