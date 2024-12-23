// sequence_aleatoire.c
#include "projet.h"

char generer_nucleotide() {
    double rand_num = (double)rand() / RAND_MAX;

    if (rand_num < 0.246) return 'A';  // A : 24.6%
    if (rand_num < 0.490) return 'T';  // T : 24.4%
    if (rand_num < 0.745) return 'C';  // C : 25.5%
    return 'G';                        // G : 25.5%
}

char* generer_sequence(int longueur_seq) {
    char *sequence = malloc((longueur_seq + 1) * sizeof(char));
    if (sequence == NULL) {
        printf("Erreur : Allocation mémoire échouée\n");
        return NULL;
    }
    for (int i = 0; i < longueur_seq; i++) {
        sequence[i] = generer_nucleotide();
    }
    sequence[longueur_seq] = '\0';
    return sequence;
}

void ecrire_fichier_sequence(const char *sequence, const char *chemin_fichier) {
    FILE *fichier = fopen(chemin_fichier, "w");
    if (fichier == NULL) {
        perror("Erreur lors de l'ouverture du fichier pour écriture");
        return;
    }
    fprintf(fichier, ">Séquence_ADN_aléatoire\n");
    int longueur = strlen(sequence);
    for (int i = 0; i < longueur; i++) {
        fputc(sequence[i], fichier);
        if ((i + 1) % 80 == 0) fputc('\n', fichier);
    }
    if (longueur % 80 != 0) fputc('\n', fichier);
    fclose(fichier);
    printf("Séquence écrite dans %s\n", chemin_fichier);
}
