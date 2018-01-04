// nc_IO.cc

#include <iostream> // For IO to command line
#include <string> // c++ strings
#include <time.h> // time_t, time, ctime
#include <netcdf.h> // For NetCDF interface
#include "nc_IO.h" // Include header file for consistency check

/*
    Module containing IO functionalities with NetCDF.
    https://www.unidata.ucar.edu/software/netcdf/netcdf-4/newdocs/netcdf-c/
*/

// Internal routine declarations
void ErrorHandler(const int status);
std::string DateTime();
// ***************
// ***************

newDS_t::newDS_t(const std::string GlobalAttr, 
    const int dims_num_, const std::string*const dim_names, const int*const dim_lengths,
    const int vars_num_, const std::string*const var_names, const std::string path)
    :dims_num(dims_num_), dimid(new int [dims_num_]), coord_varid(new int [dims_num_]), 
    vars_num(vars_num_), varid(new int [vars_num_])
{
    std::cout << "Instance of newDS_t created" << std::endl;
    /* The dataset is created; existing files of the same name are overwritten. 
    The dataset ID is assigned to ncid. Note that we don't use NetCDF4 to accomodate 
    the python interface. */
    std::string Filename = path + DateTime();
    ErrorHandler( nc_create(Filename.c_str(), NC_CLOBBER, &ncid) );
    
    // The string GlobalAttr is set as a global attribute of the dataset.
    const size_t len = GlobalAttr.length() + 1; // Length of char array (with terminator)
    ErrorHandler( nc_put_att_text(ncid, NC_GLOBAL, "GlobalAttributes", len, 
                                  GlobalAttr.c_str()) );
    
    /* The dimensions are created. Their names and lengths are supplied as arguments. */
    for (int i=0; i<dims_num; ++i)
        ErrorHandler( nc_def_dim(ncid, dim_names[i].c_str(), dim_lengths[i], 
                                 &(dimid[i])) );
    
    // Define coordinate variables
    for (int i=0; i<dims_num; ++i)
        ErrorHandler( nc_def_var(ncid, dim_names[i].c_str(), NC_FLOAT, 1, &(dimid[i]), 
                                 &(coord_varid[i])) );
    
    // Define data variables. We presume they each span all dimensions
    for (int j=0; j<vars_num; ++j)
        ErrorHandler( nc_def_var(ncid, var_names[j].c_str(), NC_FLOAT, dims_num, dimid, 
                                 &(varid[j])) );
    
    // Exit define mode and enter data mode.
    ErrorHandler( nc_enddef(ncid) );
}

void newDS_t::WriteCoordVars(const double*const*const coord_vars)
{
    // Write coordinate variables
    // The NetCDF C interface expects row-major layout.
    for (int i=0; i<dims_num; ++i)
        ErrorHandler( nc_put_var_double(ncid, coord_varid[i], coord_vars[i]) );
}

void newDS_t::WriteVars(const double*const*const vars)
{
    // Write variables
    // The NetCDF C interface expects row-major layout.
    for (int j=0; j<vars_num; ++j)
        ErrorHandler( nc_put_var_double(ncid, varid[j], vars[j]) );
}

newDS_t::~newDS_t()
{
    std::cout << "Instance of newDS_t deleted" << std::endl;
    
    /* Closes the datased referred to by ncid. */
    ErrorHandler( nc_close(ncid) );
    
    // Deallocate memory
    delete [] dimid;
    delete [] coord_varid;
    delete [] varid;
}


// ***************

// Internal routines
void ErrorHandler(const int status)
{
    /*
        INTERNAL ROUTINE
        Handles NetCDF errors using the function nc_strerror() from the NetCDF C 
        interface. If the value of retval indicates an error, a message is output.
    */
    if (status != NC_NOERR)
    {
        std::cerr << "\n***NetCDF error***\n"
                  << nc_strerror(status) << "\n\n";
    }
}

std::string DateTime()
{
    /* Returns a c++ string with the current local time and date in the specified format.
    http://www.cplusplus.com/reference/ctime/strftime/ */
    time_t rawtime;
    struct tm* timeinfo;
    char buffer [100];
    
    time(&rawtime); // Assign time to rawtime
    timeinfo = localtime(&rawtime); // Use this to assign timeinfo
    
    // Use this to assign the time with the proper format to buffer
    strftime(buffer,100,"%Y-%m-%d_%Hh%Mmin%Ss",timeinfo);
    
    // Define a c++ string from buffer
    std::string TimeString = buffer;
    
    return TimeString;
}
