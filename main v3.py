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
    all_genes = []
    all_markers = []
    rf_matrix = []

    # Add all the genes and markers to seperate lists.
    for gene, markers in dictionary.items():
        all_genes.append(gene)
        all_markers.append(markers)

    # Make a pairwise alignment of the markers of all the genes and
    # count the number of times certain combinations occur. Calculate
    # the percentages of these combinations. Add the names of the genes
    # and the percentages to a dictionary.
    for i in range(len(all_markers)):
        rf_percentages = {}
        for k in range(len(all_markers)):
            recombinations = 0
            total_count = 0
            for n in range(len(all_markers[i])):
                if str(all_markers[i][n] + all_markers[k][n]) == "ab":
                    recombinations += 1
                if str(all_markers[i][n] + all_markers[k][n]) == "ba":
                    recombinations += 1
                if "-" not in str(all_markers[i][n] +
                                  all_markers[k][n]):
                    total_count += 1
            rf_perc = round(recombinations / total_count * 100, 1)
            rf_percentages[all_genes[k]] = rf_perc
        rf_matrix.append(rf_percentages)

    # Sort the dictionaries.
    for i in range(len(rf_matrix)):
        temp_dict = sorted(rf_matrix[i].items(), key=lambda x: x[1])
        sorted_dict = dict(temp_dict)
        rf_matrix[i] = sorted_dict

    # Write the results to a file.
    with open("mapchart.txt", "w") as rf_open:
        for i in range(len(all_genes)):
            rf_open.write("Group " + all_genes[i])
            rf_open.write("\n")
            for key, value in rf_matrix[i].items():
                rf_open.write(key + ": " + str(value))
                rf_open.write("\n")
            rf_open.write("\n")

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
        plant = []
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
    calc_rf(gene_markers)
    print("*** DONE ***")
