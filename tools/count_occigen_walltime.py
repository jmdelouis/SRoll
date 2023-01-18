#!/usr/bin/env python

from __future__ import division, print_function

# sample occigen sacct line:
# $ sacct -S 2020-05-01 -X -u smottet -o start,jobid,jobname%-18,elapsed,cputime%12,reserved,state,reqTRES%-32,aveRSS,partition

sacct = """
2020-05-07T12:31:13 10264538     APR20_100-1a_000     01:42:53  80-14:12:24   00:02:04  COMPLETED cpu=512,mem=1100G,node=22                          all 
2020-05-07T14:42:34 10267068     APR20_100-1b_000     01:36:46  82-18:56:32   00:00:01  COMPLETED cpu=512,mem=1100G,node=22                          all 
2020-05-07T16:57:19 10269083     APR20_100-4a_000     01:43:41  88-16:57:52   00:00:02  COMPLETED cpu=512,mem=1100G,node=22                          all 
2020-05-07T16:57:23 10269084     APR20_100-4b_000     01:44:14  89-04:15:28   00:00:02  COMPLETED cpu=512,mem=1100G,node=22                          all 
2020-05-07T19:07:07 10272119     APR20_100-2a_000     01:56:12  99-09:58:24   00:10:35  COMPLETED cpu=512,mem=1100G,node=22                          all 
2020-05-07T19:07:07 10272120     APR20_100-2b_000     01:55:36  98-21:39:12   00:10:31  COMPLETED cpu=512,mem=1100G,node=22                          all 
2020-05-07T19:07:07 10272122     APR20_100-3a_000     01:52:56  96-14:53:52   00:10:26  COMPLETED cpu=512,mem=1100G,node=22                          all 
2020-05-07T19:08:07 10272123     APR20_100-3b_000     01:54:10  97-16:13:20   00:11:22  COMPLETED cpu=512,mem=1100G,node=22                          all 
2020-05-08T09:49:37 10274934     APR20_143-1a_000     01:43:25  88-11:29:20   00:00:00  COMPLETED cpu=512,mem=1100G,node=22                          all 
2020-05-08T09:50:03 10274945     APR20_143-1b_000     01:44:01  88-23:48:32   00:00:01  COMPLETED cpu=512,mem=1100G,node=22                          all 
2020-05-08T11:33:26 10275058     APR20_143-3a_000     01:44:21  89-06:39:12   00:00:00  COMPLETED cpu=512,mem=1100G,node=22                          all 
2020-05-08T11:33:34 10275059     APR20_143-3b_000     01:43:47  88-19:01:04   00:00:00  COMPLETED cpu=512,mem=1100G,node=22                          all 
2020-05-08T13:18:18 10275197     APR20_143-2a_000     01:41:59  87-06:03:28   00:00:01  COMPLETED cpu=512,mem=1100G,node=22                          all 
2020-05-08T13:18:23 10275198     APR20_143-2b_000     01:44:32  89-10:25:04   00:00:00  COMPLETED cpu=512,mem=1100G,node=22                          all 
2020-05-08T14:36:30 10275564     APR20_143-4a_000     01:43:29  88-12:51:28   00:00:00  COMPLETED cpu=512,mem=1100G,node=22                          all 
2020-05-08T14:36:38 10275565     APR20_143-4b_000     01:44:33  89-10:45:36   00:00:01  COMPLETED cpu=512,mem=1100G,node=22                          all 
2020-05-08T16:41:45 10275938     APR20_143-5_000      01:48:33  92-20:53:36   00:00:01  COMPLETED cpu=512,mem=1100G,node=22                          all 
2020-05-08T16:41:51 10275939     APR20_143-6_000      01:49:24  91-04:00:00   00:00:01  COMPLETED cpu=512,mem=1100G,node=22                          all 
2020-05-08T16:41:55 10275940     APR20_143-7_000      01:49:37  93-18:47:44   00:00:00  COMPLETED cpu=512,mem=1100G,node=22                          all 
2020-05-08T18:45:53 10276056     APR20_217-5a_000     01:47:08  89-06:40:00   00:00:01  COMPLETED cpu=512,mem=1100G,node=22                          all 
2020-05-08T18:45:59 10276057     APR20_217-5b_000     01:45:30  87-22:00:00   00:00:01  COMPLETED cpu=512,mem=1100G,node=22                          all 
2020-05-08T20:26:57 10276182     APR20_217-7a_000     01:43:45  86-11:00:00   00:00:00  COMPLETED cpu=512,mem=1100G,node=22                          all 
2020-05-08T20:27:03 10276183     APR20_217-7b_000     01:45:10  89-23:25:20   00:00:01  COMPLETED cpu=512,mem=1100G,node=22                          all 
2020-05-08T21:58:58 10276278     APR20_217-6a_000     01:42:20  87-13:14:40   00:00:01  COMPLETED cpu=512,mem=1100G,node=22                          all 
2020-05-08T21:59:03 10276279     APR20_217-6b_000     01:43:44  88-17:59:28   00:00:00  COMPLETED cpu=512,mem=1100G,node=22                          all 
2020-05-08T23:01:34 10276370     APR20_217-8a_000     01:44:51  89-16:55:12   00:00:01  COMPLETED cpu=512,mem=1100G,node=22                          all 
2020-05-08T23:01:39 10276371     APR20_217-8b_000     01:44:30  89-09:44:00   00:00:00  COMPLETED cpu=512,mem=1100G,node=22                          all 
2020-05-09T12:20:22 10276904     APR20_217-1_000      01:47:19  89-10:20:00   00:00:00  COMPLETED cpu=512,mem=1100G,node=22                          all 
2020-05-09T12:20:27 10276905     APR20_217-2_000      01:46:53  89-01:40:00   00:00:00  COMPLETED cpu=512,mem=1100G,node=22                          all 
2020-05-09T14:04:51 10277371     APR20_217-3_000      01:46:32  91-03:29:04   00:00:01  COMPLETED cpu=512,mem=1100G,node=22                          all 
2020-05-09T14:04:57 10277372     APR20_217-4_000      01:46:29  91-02:27:28   00:00:01  COMPLETED cpu=512,mem=1100G,node=22                          all 
2020-05-09T15:25:43 10277469     APR20_353-3a_000     01:45:00  89-20:00:00   00:00:01  COMPLETED cpu=512,mem=1100G,node=22                          all 
2020-05-09T15:25:48 10277470     APR20_353-3b_000     01:43:49  88-19:42:08   00:00:01  COMPLETED cpu=512,mem=1100G,node=22                          all 
2020-05-09T17:20:02 10277562     APR20_353-5a_000     01:44:51  89-16:55:12   00:00:02  COMPLETED cpu=512,mem=1100G,node=22                          all 
2020-05-09T17:20:06 10277563     APR20_353-5b_000     01:44:51  89-16:55:12   00:00:01  COMPLETED cpu=512,mem=1100G,node=22                          all 
2020-05-09T18:50:25 10277641     APR20_353-4a_000     01:47:49  89-20:20:00   00:00:03  COMPLETED cpu=512,mem=1100G,node=22                          all 
2020-05-09T18:50:29 10277642     APR20_353-4b_000     01:47:52  89-21:20:00   00:00:01  COMPLETED cpu=512,mem=1100G,node=22                          all 
2020-05-09T20:12:53 10277731     APR20_353-6a_000     01:45:50  88-04:40:00   00:00:00  COMPLETED cpu=512,mem=1100G,node=22                          all 
2020-05-09T20:12:58 10277732     APR20_353-6b_000     01:45:04  87-13:20:00   00:00:00  COMPLETED cpu=512,mem=1100G,node=22                          all 
2020-05-09T22:38:13 10277942     APR20_353-1_000      01:44:54  87-10:00:00   00:00:02  COMPLETED cpu=512,mem=1100G,node=22                          all 
2020-05-09T22:38:17 10277943     APR20_353-2_000      01:45:32  87-22:40:00   00:00:00  COMPLETED cpu=512,mem=1100G,node=22                          all 
2020-05-09T23:21:35 10277984     APR20_353-7_000      01:48:25  90-08:20:00   00:00:01  COMPLETED cpu=512,mem=1100G,node=22                          all 
2020-05-09T23:21:42 10277985     APR20_353-8_000      01:47:36  89-16:00:00   00:00:00  COMPLETED cpu=512,mem=1100G,node=22                          all 
2020-05-10T10:22:02 10278484     APR20_545-1_000      00:09:09   7-15:00:00   00:00:01  COMPLETED cpu=512,mem=1100G,node=22                          all 
2020-05-10T10:25:22 10278490     APR20_545-2_000      00:07:06   5-22:00:00   00:00:00  COMPLETED cpu=512,mem=1100G,node=22                          all 
2020-05-10T10:25:30 10278491     APR20_545-4_000      00:07:34   6-07:20:00   00:00:01  COMPLETED cpu=512,mem=1100G,node=22                          all 
2020-05-10T11:40:49 10278609     APR20_857-1_000      00:08:11   6-19:40:00   00:00:01  COMPLETED cpu=512,mem=1100G,node=22                          all 
2020-05-10T11:40:54 10278610     APR20_857-2_000      00:08:00   6-16:00:00   00:00:01  COMPLETED cpu=512,mem=1100G,node=22                          all 
2020-05-10T11:40:58 10278611     APR20_857-3_000      00:07:02   5-20:40:00   00:00:00  COMPLETED cpu=512,mem=1100G,node=22                          all 
2020-05-10T11:41:03 10278612     APR20_857-4_000      00:08:07   6-18:20:00   00:00:01  COMPLETED cpu=512,mem=1100G,node=22                          all 
"""
sacct = """
2020-05-17T10:42:07 10321618     S21APR20_000_100g+   00:49:43  71-06:15:12   00:00:01  COMPLETED cpu=1024,mem=5074000M,node=43                    hsw24 
2020-05-17T12:01:15 10321790     S21APR20_000_143g+   00:58:56  84-11:18:24   00:00:01  COMPLETED cpu=1024,mem=5074000M,node=43                    hsw24 
2020-05-17T12:01:26 10321791     S21APR20_000_217g+   01:20:35 115-12:04:00   00:00:01  COMPLETED cpu=1024,mem=5074000M,node=43                    hsw24 
2020-05-17T12:01:37 10321792     S21APR20_000_353g+   02:23:15 205-07:48:00   00:00:01  COMPLETED cpu=1024,mem=5074000M,node=43                    hsw24 
2020-05-17T12:01:46 10321793     S21APR20_000_545g+   00:14:53  10-21:56:48   00:00:01  COMPLETED cpu=512,mem=2596000M,node=22                     hsw24 
2020-05-17T12:01:52 10321794     S21APR20_000_857g+   00:28:17  20-17:47:12   00:00:01  COMPLETED cpu=512,mem=2596000M,node=22                     hsw24 
"""

max_wt = 0
tot_wt = 0
njobs  = 0
for job in sacct.splitlines():
  cols = job.split()
  if len( cols) > 4:
    njobs += 1
    walltime = cols[4]
    if "-" in walltime:
      days, hours = walltime.split("-")
    else:
      days = "0"
      hours = walltime
    h, m, s = hours.split(":")
    wt = int(days) * 24 + int( h) + int( m) / 60 + int( s) / 3600
    wt /= 2 # sacct cputime column shows multithreading, which is not taken into account in cines allocation
    tot_wt += wt
    max_wt = max( max_wt, wt)

print( "number of jobs: %d" % njobs)
print( "total walltime hours: %d" % tot_wt)
print( "max job walltime: %d" % max_wt)
print( "average walltime: %d" % (tot_wt / njobs))
