#!usr/bin/python
import math
import csv
import numpy as np
import sys


# This function is from Hill & Weir 2011 'Variation in actual Relationship as a consequence of Mendelian sampling and
# linkage'. The formulae to calculate the estimated variation in the proportion of genome shared by decent in a
# relationship are all in terms of this function

def phi(n, l):

    if n == 0:
        return 0

    else:
        y = 0
        for r in range(1, n+1):

            y += (math.factorial(n)/(math.factorial(r) * math.factorial(n-r))) * (((2 * r * l) - 1 + math.exp(-2 * r * l)) / (r ** 2))

        return (1. / (2 * (l ** 2) * (4 ** n))) * y


# These functions are the formulae for the variation in the relationship coefficient for the different types of
# relationships. The l is the map length of the chromosomes in question. g = log2(R), this is the generations between
# the individuals (see Hill & Weir 2011)

def lineal(g, l):
    return phi(g - 1, l) / 4


def half_sib(g, l):
    return (4 * phi(g, l) - 2 * phi(g - 1, l) + (phi(g - 2, l) / 2)) / 4


def uncle_nephew(g, l):
    return (8 * phi(g + 1, l) - 4 * phi(g, l) + (2 * phi(g - 1, l) + phi(g - 2, l)) / 4) / 4


def cousins(g, l):
    return (8 * phi(g + 1, l) - 4 * phi(g, l) + (3 / 2) * phi(g - 1, l) - phi(g - 2, l) / 2 + phi(g - 3, l) / 8) / 4


def full_siblings(l):
    return 2 * phi(2, l) - phi(1, l)


# double first cousins
def dfc(l):
    return 4 * phi(4, l) - 2 * phi(3, l) + (3 / 4) * phi(2, l) - (1 / 4) * phi(1, l)


# This calculates the ibd matrix from a fam file, with individuals ordered by decreasing age

def ibd():

    # open file specified in command line argument
    with open(sys.argv[1], 'r') as rfile:
        r = csv.reader(rfile)

        # In the script I use the individuals position as it's id, therefore I have to save the original id with it's
        # corresponding number so that I can return the original ids with the result
        orig_ind = {}
        ind = {}
        count = 1

        # I use '0' as a missing individual, I'll have to specify the format the input file will have to be
        orig_ind['0'] = 0

        # Here I save each individual's id and their parent's ids, I first convert the ids into a numeric orderfirst so
        # that it is easier to work with in the code
        for row in r:
            orig_ind[row[0]] = count
            ind[count] = [orig_ind[row[1]], orig_ind[row[2]]]
            count += 1

        # initialise the grm with 0's

        grm = np.array([[0.0]*len(ind)] * len(ind))

        # I use an algorithm to calculate the grm, this can be found online and in a lot of books ( I can't remember
        # the names)
        for i in sorted(ind.iterkeys()):
            for j in range(0, i - 1):
                if ind[i][0] != 0 and ind[i][1] != 0:
                    grm[i-1][j] = (grm[ind[i][0] - 1][j] + grm[ind[i][1] - 1][j]) / 2

                else:
                    if ind[i][0] != 0 and ind[i][1] == 0:
                        grm[i - 1][j] = (grm[ind[i][0] - 1][j]) / 2
                    elif ind[i][0] == 0 and ind[i][1] != 0:
                        grm[i - 1][j] = (grm[ind[i][1] - 1][j]) / 2
                    else:
                        grm[i - 1][j] = 0

                grm[j][i - 1] = grm[i - 1][j]

            if ind[i][0] != 0 and ind[i][1] != 0:
                grm[i-1][i-1] = 1 + grm[ind[i][0] - 1][ind[i][1] - 1] / 2

            else:
                grm[i - 1][i - 1] = 1

    # Initialise the variance matrix and a matrix storing the relationships
    var = np.array([[0.0]*len(ind)] * len(ind))

    rel = np.array([[''] * len(ind)] * len(ind))

    # Looping over the individuals
    for a in range(0, len(grm[0, ])):
        # Then loop over the individuals before the individual in question in the pedigree (i.e. older than this
        # individual)
        for b in range(0, a):

            # If the R = 1/2, The individuals are either parent-offspring or full siblings

            if grm[a][b] == 0.5:
                if b+1 in ind[a+1]:
                    var[a][b] = lineal(1, 1.5)
                    rel[a][b] = 'P'

                else:
                    var[a][b] = full_siblings(1.5)
                    rel[a][b] = 'F'
                rel[b][a] = rel[a][b]

            # If the R = 1/4, the individuals are either half siblings, Uncle-Nephew, double first cousins, or
            # grandparent grandoffspring
            elif grm[a][b] == 0.25:

                if set(ind[a+1]).intersection(set(ind[b+1])) != set() and set(ind[a+1]).intersection(set(ind[b+1])) != set([0]):
                    var[a][b] = half_sib(2, 1.5)
                    rel[a][b] = 'H'

                elif rel[ind[a+1][0]][b+1] == 'F' or rel[ind[a+1][1]][b+1] == 'F':
                    var[a][b] = uncle_nephew(2, 1.5)
                    rel[a][b] = 'U'

                elif (rel[ind[a+1][0]][ind[b+1][0]] == 'F' and rel[ind[a+1][1]][ind[b+1][1]] == 'F') or (rel[ind[a+1][1]][ind[b+1][0]] == 'F' and rel[ind[a+1][0]][ind[b+1][1]] == 'F'):
                    var[a][b] = dfc(1.5)
                    rel[a][b] = 'D'

                else:
                    var[a][b] = lineal(2, 1.5)
                    rel[a][b] = 'G'

                rel[b][a] = rel[a][b]

            # If R=0 they are not related
            elif grm[a][b] == 0:

                var[a][b] = 0

            # If R < 1/4, then the individuals are descendants of one of the above relationships
            else:
                # We want to check the younger individuals ancestors for a relationship recorded above to the older
                # individual
                trip = 0
                check = []

                # We want to make sure we are not check for missing individuals...
                if ind[a + 1][0] != 0:
                    check.append(ind[a + 1][0])

                if ind[a + 1][1] != 0:
                    check.append(ind[a + 1][1])

                k = 0

                while trip == 0:

                    # If we find a relationship we record it and exit the loop
                    if rel[check[k] - 1][b] != '':
                        rel[a][b] = rel[check[k] - 1][b]
                        rel[b][a] = rel[a][b]
                        trip = 1

                    # else we add the parents of the individual being checked against the original to the queue of
                    # people to be checked
                    else:
                        if ind[check[k]][0] != 0:
                            check.append(ind[check[k]][0])

                        if ind[check[k]][1] != 0:
                            check.append(ind[check[k]][1])
                    k += 1

                # once we find a relationship we can calculate the variance since we know that g = log2(R)
                g = int(abs(math.log(grm[a][b], 2)))

                if rel[a][b] == 'H':
                    var[a][b] = half_sib(g, 1.5)

                elif rel[a][b] == 'U':
                    var[a][b] = uncle_nephew(g, 1.5)

                elif rel[a][b] == 'G':
                    var[a][b] = lineal(g, 1.5)

                elif rel[a][b] == 'F':
                    var[a][b] = cousins(g, 1.5)

                var[b][a] = var[a][b]

    # We make sure to return the original ids instead of the numerical order
    inv_orig_ind = {v: k for k, v in orig_ind.iteritems()}

    return [grm, var, inv_orig_ind]

result = ibd()

f = open('temp.txt', 'w')
x = csv.writer(f, delimiter='\t')

x.writerow(['individual_1', 'individual_2', 'Relationship_coefficient', 'SD'])

for i in range(0, len(result[0])):
    for j in range(0, i):
        if result[1][i][j] != 0:
            x.writerow([result[2][i + 1], result[2][j + 1], result[0][i][j], result[1][i][j] / 2])

f.close()