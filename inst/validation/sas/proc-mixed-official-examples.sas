/*
  Reproducible PROC MIXED specifications used by lmmix validation.

  Source: SAS Institute Inc. (2017), SAS/STAT 14.3 User's Guide,
  Chapter 79, Examples 79.1, 79.2, and 79.5.
  https://documentation.sas.com/api/docsets/statug/14.3/content/mixed.pdf?locale=en
*/

data split_plot;
  input Block A B Y @@;
  datalines;
1 1 1 56  1 1 2 41  1 2 1 50  1 2 2 36  1 3 1 39  1 3 2 35
2 1 1 30  2 1 2 25  2 2 1 36  2 2 2 28  2 3 1 33  2 3 2 30
3 1 1 32  3 1 2 24  3 2 1 31  3 2 2 27  3 3 1 15  3 3 2 19
4 1 1 30  4 1 2 25  4 2 1 35  4 2 2 30  4 3 1 17  4 3 2 18
;

ods output CovParms=split_cov FitStatistics=split_fit Tests3=split_tests;
proc mixed data=split_plot;
  class A B Block;
  model Y = A B A*B;
  random Block A*Block;
run;

data growth;
  input Person Gender $ y1 y2 y3 y4;
  y=y1; Age=8; output;
  y=y2; Age=10; output;
  y=y3; Age=12; output;
  y=y4; Age=14; output;
  drop y1-y4;
  datalines;
 1 F 21.0 20.0 21.5 23.0
 2 F 21.0 21.5 24.0 25.5
 3 F 20.5 24.0 24.5 26.0
 4 F 23.5 24.5 25.0 26.5
 5 F 21.5 23.0 22.5 23.5
 6 F 20.0 21.0 21.0 22.5
 7 F 21.5 22.5 23.0 25.0
 8 F 23.0 23.0 23.5 24.0
 9 F 20.0 21.0 22.0 21.5
10 F 16.5 19.0 19.0 19.5
11 F 24.5 25.0 28.0 28.0
12 M 26.0 25.0 29.0 31.0
13 M 21.5 22.5 23.0 26.5
14 M 23.0 22.5 24.0 27.5
15 M 25.5 27.5 26.5 27.0
16 M 20.0 23.5 22.5 26.0
17 M 24.5 25.5 27.0 28.5
18 M 22.0 22.0 24.5 26.5
19 M 24.0 21.5 24.5 25.5
20 M 23.0 20.5 31.0 26.0
21 M 27.5 28.0 31.0 31.5
22 M 23.0 23.0 23.5 25.0
23 M 21.5 23.5 24.0 28.0
24 M 17.0 24.5 26.0 29.5
25 M 22.5 25.5 25.5 26.0
26 M 23.0 24.5 26.0 30.0
27 M 22.0 21.5 23.5 25.0
;

ods output CovParms=growth_un_cov FitStatistics=growth_un_fit
  SolutionF=growth_un_fixed Tests3=growth_un_tests;
proc mixed data=growth method=ml covtest;
  class Person Gender;
  model y = Gender Age Gender*Age / s;
  repeated / type=un subject=Person r;
run;

ods output CovParms=growth_cs_cov FitStatistics=growth_cs_fit
  SolutionF=growth_cs_fixed Tests3=growth_cs_tests;
proc mixed data=growth method=ml covtest;
  class Person Gender;
  model y = Gender Age Gender*Age / s;
  repeated / type=cs subject=Person r;
run;

data random_coefficients;
  input Batch Month @@;
  do Replicate = 1 to 6;
    input Y @@;
    output;
  end;
  datalines;
1 0 101.2 103.3 103.3 102.1 104.4 102.4
1 1 98.8 99.4 99.7 99.5 . .
1 3 98.4 99.0 97.3 99.8 . .
1 6 101.5 100.2 101.7 102.7 . .
1 9 96.3 97.2 97.2 96.3 . .
1 12 97.3 97.9 96.8 97.7 97.7 96.7
2 0 102.6 102.7 102.4 102.1 102.9 102.6
2 1 99.1 99.0 99.9 100.6 . .
2 3 105.7 103.3 103.4 104.0 . .
2 6 101.3 101.5 100.9 101.4 . .
2 9 94.1 96.5 97.2 95.6 . .
2 12 93.1 92.8 95.4 92.2 92.2 93.0
3 0 105.1 103.9 106.1 104.1 103.7 104.6
3 1 102.2 102.0 100.8 99.8 . .
3 3 101.2 101.8 100.8 102.6 . .
3 6 101.1 102.0 100.1 100.2 . .
3 9 100.9 99.5 102.2 100.8 . .
3 12 97.8 98.3 96.9 98.4 96.9 96.5
;

ods output CovParms=random_cov FitStatistics=random_fit
  SolutionF=random_fixed SolutionR=random_effects Tests3=random_tests;
proc mixed data=random_coefficients;
  class Batch;
  model Y = Month / s;
  random Int Month / type=un subject=Batch s;
run;
