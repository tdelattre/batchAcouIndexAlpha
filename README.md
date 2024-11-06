# batchAcouIndexAlpha
batch computing acoustic indexes using parallel processing

# Objectif du script : 
Appliquer les indices acoustiques sur les fichiers wav bruts d'une minutes collectés par des capteurs wildlife acoustics, sur la base d'une structure de fichiers

``` projet / session / parcelle / Data2 / fichiers sons ```

Le script consiste essentiellement en une sérialisation de la fonction AcouIndexAlphaCustom créée par Amandine Gasc, qui calcule les indices acoustiques les plus communs (ACI, NDSI, BI, ADI, NP).
Il est conçu pour traiter des fichiers wav d'une minute, mais cela est facilement adaptable.

Cette fonction peut être retrouvée ici : https://github.com/agasc/Soundscape-analysis-with-R/blob/master/AcouIndexAlpha.r
Cette fonction est dans le fichier AcouIndexAlphaCustomSilent.r, merci de télécharger les deux fichiers pour que l'ensemble fonctionne.

Le script tire partie des bibliothèques de parallel processing de R et devrait par conséquent être raisonnablement portable. Cependant ces techniques reposent fortement sur des spécificités matérielles qui peuvent varier d'une machine à l'autre. Merci de me signaler tout discfonctionnement pour aider à améliorer le code et le rendre disponible au plus grand nombre.

