# do some analysis of variance
import sys
import csv
from tabulate import tabulate
from scipy import stats

ANOVA_HEADERS = ['src', 'DF', 'SS', 'MS', 'F', 'p']


def sum_x(data):
    total = 0
    for row in data:
        total += sum(row)
    return total

def sum_xsq(data):
    total = 0
    for row in data:
        for item in row:
            total += item**2
    return total


def prep_data(file, reps):
    reader = csv.reader(file)
    init_data = []
    for line in reader:
        init_data.append([float(item) for item in line])

    n_rows = len(init_data)
    n_cols = len(init_data[0])
    n = n_cols * n_rows

    gridsum_data = []
    for i in range(0, n_rows, reps):
        summed_line = []
        for j in range(n_cols):
            block_sum = 0
            for k in range(i, i+reps, 1):
                block_sum += init_data[k][j]
            summed_line.append(block_sum)
        gridsum_data.append(summed_line)

    col_sums = [0 for _ in range(n_cols)]
    row_sums = []
    for row in gridsum_data:
        row_sums.append(sum(row))
        for i in range(len(row)):
            col_sums[i] += row[i]

    grand_total = sum(col_sums)

    return init_data, gridsum_data, col_sums, row_sums, grand_total, n, n_cols, n_rows


def pretty_anova_tbl(anova, sources):
    anova_tbl = []
    for src in sources:
        line = [src, ]
        for attr in ANOVA_HEADERS[1:]:
            if attr in anova[src]:
                line.append(anova[src][attr])
            else:
                line.append(None)
        anova_tbl.append(line)

    print(tabulate(anova_tbl, headers=ANOVA_HEADERS))


def simple_anova(file, reps):

    data, _, treat_sums, _, grand_total, n, treat_n, obs_per_treat = prep_data(file, reps)

    anova = {}; sources = ['treats', 'err', 'total']
    for item in sources:
        anova[item] = {}

    # DF's
    anova['treats']['DF'] = treat_n-1
    anova['total']['DF'] = n-1
    anova['err']['DF'] = anova['total']['DF'] - anova['treats']['DF']

    # SS's
    anova['total']['SS'] = sum_xsq(data) - ((sum_x(data)**2)/n)
    ssTreats = 0
    for i in range(treat_n):
        ssTreats += (treat_sums[i]**2)/obs_per_treat
    ssTreats -= (grand_total**2)/n
    anova['treats']['SS'] = ssTreats
    anova['err']['SS'] = anova['total']['SS'] - anova['treats']['SS']

    # MS's
    anova['treats']['MS'] = anova['treats']['SS']/anova['treats']['DF']
    anova['err']['MS'] = anova['err']['SS']/anova['err']['DF']

    # F
    anova['treats']['F'] = anova['treats']['MS']/anova['err']['MS']

    # p
    anova['treats']['p'] = stats.f(anova['treats']['DF'], anova['err']['DF']).sf(anova['treats']['F'])

    pretty_anova_tbl(anova, sources)


def blocked_anova(file, reps):

    data, _, treat_sums, block_sums, grand_total, n, n_treats, obs_per_treat = prep_data(file, reps)
    n_blocks = obs_per_treat//reps
    obs_per_block = n // n_blocks

    anova = {}; sources = ['treats', 'blocks', 'err', 'total']
    for item in sources:
        anova[item] = {}

    # DF's
    anova['total']['DF'] = n-1
    anova['treats']['DF'] = n_treats-1
    anova['blocks']['DF'] = n_blocks-1
    anova['err']['DF'] = anova['total']['DF'] - anova['treats']['DF'] - anova['blocks']['DF']

    # SS's
    anova['total']['SS'] = sum_xsq(data) - ((sum_x(data) ** 2) / n)
    ssTreats = 0
    for i in range(n_treats):
        ssTreats += (treat_sums[i] ** 2) / obs_per_treat
    ssTreats -= (grand_total ** 2) / n
    anova['treats']['SS'] = ssTreats
    ssBlocks = 0
    for i in range(n_blocks):
        ssBlocks += (block_sums[i]**2) / obs_per_block
    ssBlocks -= (grand_total**2) / n
    anova['blocks']['SS'] = ssBlocks
    anova['err']['SS'] = anova['total']['SS'] - anova['treats']['SS'] - anova['blocks']['SS']

    # MS's
    anova['treats']['MS'] = anova['treats']['SS'] / anova['treats']['DF']
    anova['blocks']['MS'] = anova['blocks']['SS'] / anova['blocks']['DF']
    anova['err']['MS'] = anova['err']['SS'] / anova['err']['DF']

    # F
    anova['treats']['F'] = anova['treats']['MS'] / anova['err']['MS']
    anova['blocks']['F'] = anova['blocks']['MS'] / anova['err']['MS']

    # p
    anova['treats']['p'] = stats.f(anova['treats']['DF'], anova['err']['DF']).sf(anova['treats']['F'])
    anova['blocks']['p'] = stats.f(anova['blocks']['DF'], anova['err']['DF']).sf(anova['blocks']['F'])

    pretty_anova_tbl(anova, sources)


def twoway_anova(file, reps):
    data, gridsum_data, fac1_sums, fac2_sums, grand_total, n, n_fac1s, obs_per_fac1 = prep_data(file, reps)
    n_fac2s = obs_per_fac1 // reps
    obs_per_fac2 = n // n_fac2s
    obs_per_square = obs_per_fac1 // n_fac2s

    anova = {}; sources = ['faccol', 'facrow', 'inter', 'err', 'total']
    for item in sources:
        anova[item] = {}

    # DF's
    anova['total']['DF'] = n - 1
    anova['faccol']['DF'] = n_fac1s-1
    anova['facrow']['DF'] = n_fac2s-1
    anova['inter']['DF'] = anova['faccol']['DF'] * anova['facrow']['DF']
    anova['err']['DF'] = anova['total']['DF'] - anova['faccol']['DF'] - anova['facrow']['DF'] - anova['inter']['DF']

    # SS's
    anova['total']['SS'] = sum_xsq(data) - ((sum_x(data) ** 2) / n)
    SSFac1 = 0
    for sum in fac1_sums:
        SSFac1 += (sum**2)/obs_per_fac1
    SSFac1 -= grand_total**2 / n
    anova['faccol']['SS'] = SSFac1
    SSFac2 = 0
    for sum in fac2_sums:
        SSFac2 += (sum**2)/obs_per_fac2
    SSFac2 -= grand_total**2 / n
    anova['facrow']['SS'] = SSFac2
    SSInter = 0
    for i in range(n_fac2s):  # rows
        for j in range(n_fac1s):  # cols
            SSInter += gridsum_data[i][j]**2 / obs_per_square
    SSInter -= grand_total**2 / n
    SSInter -= SSFac1 + SSFac2
    anova['inter']['SS'] = SSInter
    anova['err']['SS'] = anova['total']['SS'] - anova['faccol']['SS'] - anova['facrow']['SS'] - anova['inter']['SS']

    # MS's
    anova['faccol']['MS'] = anova['faccol']['SS'] / anova['faccol']['DF']
    anova['facrow']['MS'] = anova['facrow']['SS'] / anova['facrow']['DF']
    anova['inter']['MS'] = anova['inter']['SS'] / anova['inter']['DF']
    anova['err']['MS'] = anova['err']['SS'] / anova['err']['DF']

    # F
    anova['faccol']['F'] = anova['faccol']['MS'] / anova['err']['MS']
    anova['facrow']['F'] = anova['facrow']['MS'] / anova['err']['MS']
    anova['inter']['F'] = anova['inter']['MS'] / anova['err']['MS']

    # p
    anova['faccol']['p'] = stats.f(anova['faccol']['DF'], anova['err']['DF']).sf(anova['faccol']['F'])
    anova['facrow']['p'] = stats.f(anova['facrow']['DF'], anova['err']['DF']).sf(anova['facrow']['F'])
    anova['inter']['p'] = stats.f(anova['inter']['DF'], anova['err']['DF']).sf(anova['inter']['F'])

    pretty_anova_tbl(anova, sources)


if __name__ == '__main__':
    if len(sys.argv) < 4:
        print('Usage: python3', sys.argv[0], '[simple|blocked|twoway]', '[csv_file]', '[reps]')
        exit()

    method = sys.argv[1]
    reps = int(sys.argv[3])
    with open(sys.argv[2]) as file:
        if method == 'simple':
            simple_anova(file, reps)
        elif method == 'blocked':
            blocked_anova(file, reps)
        elif method == 'twoway':
            twoway_anova(file, reps)
        else:
            print('Unknown method:', method)
            print('Allowed methods: simple, blocked, twoway')


