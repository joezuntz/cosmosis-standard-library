from cosmosis.datablock import option_section, names
import consistency

def setup(options):
    verbose = options.get_bool(option_section, "verbose", default=False)
    cons = consistency.cosmology_consistency(verbose)
    return cons

def execute(block, config):
    cons = config
    cosmo = names.cosmological_parameters
    
    #Create dict of all parameters that we have already
    known_parameters = {}
    for param in cons.parameters:
        if block.has_value(cosmo, param):
            known_parameters[param] = block[cosmo,param]

    if cons.verbose:
        print "Consistency relation input parameters:", ', '.join(known_parameters.keys())

    #Run the consistency checker/parameter filler-inner.
    try:
        filled_parameters = cons(known_parameters)
    #There are two possible simpler errors:
    #too many/inconsistent parameters, or not enough.
    except consistency.UnderSpecifiedModel as error:
        print "You did not set enough cosmological parameters"
        print "to be able to deduce the rest of them:"
        print error
        return 1
    except consistency.OverSpecifiedModel as error:
        print "You set inconsistent cosmological parameters:"
        print error
        return 2

    #Set or replace the new values
    for param, value in filled_parameters.items():
        block[cosmo,param] = value

    return 0


            
