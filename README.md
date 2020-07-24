# ANOVA from scratch

This repo contains code that does statistical variance analysis from good ol'
first principles. 

## Usage outline

```sh
python3 anova.py [simple|blocked|twoway] [csv_file] [reps]
```

where 

- The first argument is one of the three options specifying what kind of ANOVA 
to perform. This will essentially determine the number of rows in the final 
ANOVA table. See below
- `[csv_file]` is the raw data csv file. No headers for rows or columns. All rows and columns should have the same lengths, respectively. Formatting below.
- `[reps]` is the number of repetitions the experiment underwent. Place multiple reps
directly underneath each other in the same block. So if 1 rep has 3 rows and 3 cols,
then 2 reps should have 6 rows and 3 cols. The program will then interpret
the csv data correctly. 

## Usage detail

### Simple

Compare any number of samples. The CSV file should have a column for each sample.

### Blocked

Compare any number of samples, taking blocks into account. CSV should have a col
for each sample, and a row for each blocking level. If 2 reps for instance, then the
number of rows should be `2*blocks` (2 reps for each blocking factor).

### Two-way

Compare any number of samples, taking a second factor and their interaction into account.
Just like with blocked, the CSV should have a col for each sample, and a row for each
secondary factor level. Scales up with reps just like the blocked design. 

## Output Example   

```sh
python3 anova.py twoway emulsion.csv 2
```

    src       DF          SS          MS          F           p
    ------  ----  ----------  ----------  ---------  ----------
    faccol     2  0.631667    0.315833    15.7917    0.00406881
    facrow     1  0.00333333  0.00333333   0.166667  0.697261
    inter      2  0.0116667   0.00583333   0.291667  0.757035
    err        6  0.12        0.02
    total     11  0.766667

## Dependencies

- `tabulate` for tabulation
- `scipy` just for getting the p-val from the F-stat


## Why did I bother?

I did a sophomoric statistics course in university where we covered ANOVA. I noticed that 
it's quite a lot of repetitive work, and it seemed easily programmable. And so I 
programmed it from first principles (no fancy `R` or `pandas`) to prove to
myself I could. It also helped a great deal in tests and exams (which were 
online in 2020 thanks to COVID). Since I wrote it myself from first princips I don't 
consider it cheating. 
