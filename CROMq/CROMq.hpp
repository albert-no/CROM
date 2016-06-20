/*
   CROMq.hpp
*/
#ifndef CROMq_H
#define CROMq_H

/* 
*/

class CROM_q
{
    /* read file and split rate and run CROM

    Parameters
    ----------
    fname :: name of the file contains sequences
    R :: overall rate
    x :: number of sequences, default = 101

    */
    std::string fname;
    double R;
    int num_x;

    compute_sub_rates();

public:
    // Constructor
    CROM_q(std::string fname_in, num_x_in);

    // Destructor
    ~CROM_q();
}

#endif
