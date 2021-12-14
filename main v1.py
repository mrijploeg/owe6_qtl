def import_data(locfile_name):
    """ Imports the data from a loc-file and saves it as a dictionary.

    :param locfile_name: string - name of the loc-file
    :return dictionary: dictionary - contains a gene name as the key
    and markers in a list as the value
    """
    # Define the variables.
    counter = 0
    gene = ""
    markers = []
    temp_dictionary = {}

    # Open the file. Loop through the lines in the file. Because the
    # first few lines contain annotation, they are skipped using a
    # counter. After the annotation lines, start saving the data. First
    # is a gene and the lines after that contain markers. Save the name
    # of the gene as a string and the markers in a list. Bundle all of
    # the data in a dictionary like PVV4: [a,a,b,...].
    with open(locfile_name, "r") as locfile_open:
        for line in locfile_open:
            if counter > 6:
                if not line.startswith(" "):
                    if not temp_dictionary and not markers:
                        gene = line.split(" ")[0]
                    else:
                        temp_dictionary[gene] = markers
                        gene = line.split(" ")[0]
                        markers = []
                else:
                    for char in line:
                        if char == "a" or char == "b" or char == "-":
                            markers.append(char)
            counter += 1
        temp_dictionary[gene] = markers

    # Delete plants that have unknown markers.
    # dictionary = delete_unknown(temp_dictionary)
    dictionary = temp_dictionary

    # Return the dictionary.
    return dictionary


def delete_unknown(dictionary):
    """ Delete plants with unknown values from the list with markers.

    :param dictionary: dictionary - contains a gene name as the key
    and markers in a list as the value
    :return: dictionary: dictionary - contains a gene name as the
    key and markers in a list as the value, but plants with missing
    values are removed
    """
    # Define the variables.
    missing = []

    # Loop through the markers in the dictionary. Save the positions
    # of the plants that are missing markers. Missing markers are
    # indicated by the "-"-symbol.
    for markers in dictionary.values():
        for marker in markers:
            if marker == "-" and markers.index(marker) not in missing:
                missing.append(markers.index(marker))

    # Print the number of missing markers.
    print("!!! MISSING MARKERS REMOVED: " + str(len(missing)) + " !!!")

    # Sort the positions of the missing markers in reversed order.
    missing.sort(reverse=True)

    # Remove the plants with any missing markers.
    for markers in dictionary.values():
        for index in missing:
            markers.pop(index)

    # Return the new dictionary.
    return dictionary


def calc_chi(dictionary):
    """ Do a chi-square test to check if the markers can be used.

    :param dictionary: dictionary - contains a gene name as the key
    and markers in a list as the value
    :return:
    """
    # Define the variables.
    gene_chi = {}
    p_value = 3.841
    delete = []

    # Loop through the genes and their markers. Calculate chi-squared
    # by comparing the real values with their expected values. If
    # chi-squared is smaller than the p-value, the real values can
    # be trusted. In this case, a = 0.05 and df = 1 are used. If the
    # chi-squared value is too high, add the gene to the list of genes
    # whose markers cannot be trusted and should be deleted.
    for gene, markers in dictionary.items():
        expected = (len(markers) / 2) - markers.count("-")
        a = markers.count("a")
        b = markers.count("b")
        chi_squared = round(((a - expected) ** 2 / expected) +
                            ((b - expected) ** 2 / expected), 3)
        gene_chi[gene] = chi_squared
        if p_value < chi_squared:
            delete.append(gene)
            print("!!! DELETED GENE: " + gene + " BECAUSE P < " +
                  str(chi_squared) + " !!!")

    # Write the results to a file.
    with open("chi_results.txt", "w") as chi_open:
        for gene, chi in gene_chi.items():
            if gene not in delete:
                chi_open.write(gene + ": " + str(chi) + "\n")
            else:
                chi_open.write(gene + ": " + str(chi) + " *\n")

    # Delete the genes that have markers that cannot be trusted.
    for gene in delete:
        dictionary.pop(gene)

    # Return the updated dictionary.
    return dictionary


def calc_rf(dictionary):
    """ Calculate the pairwise recombination frequencies and add
    those to a dictionary.

    :param dictionary: dictionary - contains a gene name as the key
    and markers in a list as the value
    :return rf_dictionary: dictionary - contains two gene names as the
    key and recombination frequencies in a list as the value
    """
    # Define the variables.
    rf_dictionary = {}
    all_genes = []
    all_markers = []

    # Add all the genes and markers to seperate lists.
    for gene, markers in dictionary.items():
        all_genes.append(gene)
        all_markers.append(markers)

    # Make a pairwise alignment of the markers of all the genes and
    # count the number of times certain combinations occur. Calculate
    # the percentages of these combinations. Add the names of the genes
    # and the percentages to a dictionary.
    for i in range(len(all_markers)):
        for k in range(len(all_markers) - 1):
            if k < i:
                continue
            else:
                aa_count = 0
                ab_count = 0
                bb_count = 0
                for n in range(len(all_markers[i])):
                    if all_markers[i][n] == "a" and \
                            all_markers[k+1][n] == "a":
                        aa_count += 1
                    elif all_markers[i][n] == "a" and \
                            all_markers[k+1][n] == "b":
                        ab_count += 1
                    elif all_markers[i][n] == "b" and \
                            all_markers[k+1][n] == "a":
                        ab_count += 1
                    elif all_markers[i][n] == "b" and \
                            all_markers[k+1][n] == "b":
                        bb_count += 1
                total_count = aa_count + ab_count + bb_count
                aa_perc = round(aa_count / total_count * 100, 1)
                ab_perc = round(ab_count / total_count * 100, 1)
                bb_perc = round(bb_count / total_count * 100, 1)
                if all_genes[i] + " + " + all_genes[k+1] \
                        not in rf_dictionary:
                    rf_dictionary[all_genes[i] + " + " +
                                  all_genes[k+1]] = \
                        [aa_perc, ab_perc, bb_perc]

    # Write the results to a file.
    with open("rf_pairwise.txt", "w") as rf_open:
        rf_open.write("GENE + GENE: [aa%, ab%, bb%]\n\n")
        for genes, freqs in rf_dictionary.items():
            rf_open.write(genes + ": " + str(freqs) + "\n")

    # Return the dictionary.
    return rf_dictionary


def create_csv(dictionary):
    """ Creates a CSV-file with all the markers per plant and the
    genes the markers correspond to.

    :param dictionary: dictionary - contains a gene name as the key
    and markers in a list as the value
    """
    # Define the variables.
    plants = {}
    genes = ["plant"]
    all_markers = []

    # Add all the genes and markers to seperate lists.
    for gene, markers in dictionary.items():
        genes.append(gene)
        all_markers.append(markers)

    # Get the markers per plant, creating lists which can later
    # be used to create a table in CSV-format.
    for i in range(162):
        plant = [str(i + 1)]
        for k in range(len(all_markers)):
            plant.append(all_markers[k][i])
        if "-" not in plant:
            plants[i] = plant
    
    # Write the data to a CSV-file.
    with open("plants.txt", "w") as plants_open:
        plants_open.write(",".join(genes))
        plants_open.write("\n")
        for key, value in plants.items():
            plants_open.write(",".join(value))
            plants_open.write("\n")


if __name__ == '__main__':
    # Define the variables.
    file_name = "CvixLer-MarkerSubset-LG1.txt"

    # Call the functions.
    print("*** IMPORTING DATA ***")
    gene_markers = import_data(file_name)
    print("*** IMPORTED DATA: " + str(len(gene_markers)) + " GENES ***")

    print("*** CHI-SQUARE TEST IN PROGRESS ***")
    gene_markers = calc_chi(gene_markers)
    print("*** GENES LEFT: " + str(len(gene_markers)) + " GENES ***")

    print("*** CREATING CSV-FILE ***")
    create_csv(gene_markers)

    print("*** CALCULATING RECOMBINATION FREQUENCIES ***")
    recombination_freqs = calc_rf(gene_markers)
    print("*** GENE PAIRS CREATED: " + str(len(recombination_freqs)) +
          " PAIRS ***")
