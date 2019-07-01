#ifndef __FEATURES_H
#define __FEATURES_H

void 
offCopy(Region *reg, BonusScores initScores);

void 
initCopy(Region *reg, Sequence seq, BonusScores initScores);

void
sdScan (Region *reg, Sequence seq, char **matrix, 
        BonusScores sdScores);

void
promoterScan (Region *reg, Sequence seq, char **matrix, 
        BonusScores promScores);

#endif /* __FEATURES_H */




