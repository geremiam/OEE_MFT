// nc_IO.h
/* Module containing IO functionalities with NetCDF. */
#ifndef NC_IO_H
#define NC_IO_H

class newDS_t
{
private:
    // Private copy constructor (prohibits copy creation)
    newDS_t(const newDS_t&);
    // Private assignment operator (prohibits assignment)
    const newDS_t& operator=(const newDS_t&);
    
    const int dims_num;
    const int vars_num;
    
    int ncid=-99;
    int*const dimid;
    int*const coord_varid;
    int*const varid;
public:
    // Constructor declaration
    newDS_t(const std::string GlobalAttr, 
      const int dims_num_, const std::string*const dim_names, const int*const dim_lengths,
      const int vars_num_, const std::string*const var_names, const std::string path="");
    // Destructor declaration
    ~newDS_t();
    
    // Write coordinate variables
    void WriteCoordVars(const double*const*const coord_vars);
    // Write variables
    void WriteVars(const double*const*const vars);
};

#endif
