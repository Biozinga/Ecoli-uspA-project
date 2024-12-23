// recherche_gene.c
#include "projet.h"

float calculer_identite(const char* seq1, const char *seq2, int longueur) {
    int matches = 0;
    for (int i = 0; i < longueur; i++) {
        if (seq1[i] == seq2[i]) {
            matches++;
        }
    }
    return (float)matches / longueur;
}
// Retourne la première position du gène
int rechercher_gene(const char* sequence_genome, const char* sequence_gene, double identite_min) {
    int longueur_genome = strlen(sequence_genome);
    int longueur_gene = strlen(sequence_gene);
    int nombre_occurence = 0;
    int position_premiere_occurence = 0;
    bool gene_trouve = false;

    for (int i = 0; i <= longueur_genome - longueur_gene; i++) {
        const char* sous_sequence = sequence_genome + i;
        float identite = calculer_identite(sous_sequence, sequence_gene, longueur_gene);

        if (identite >= identite_min) {
            printf("Gène trouvé à la position %d avec une identité de %.2f%%\n", i + 1, identite * 100);
            nombre_occurence++;
            if (!gene_trouve) {
                position_premiere_occurence = i + 1; // Position 1-indexée
                gene_trouve = true;
            }
            if (nombre_occurence > 1) {
                printf("Gène trouvé à plusieurs reprises dans le génome !\n");
                return -1;  // Sortir immédiatement si plusieurs occurrences
            }
        }
    }

    if (!gene_trouve) {
        printf("Gène non trouvé\n");
        return -1;
    }

    return position_premiere_occurence;
}
