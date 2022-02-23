#!/Users/cyrus/miniconda3/bin/python3

def read_xml(   f_name,
                fields):
    retDict = {}
    with open(f_name,'r') as f:
        text = f.read().split("\n")
    for line in text:
        for field in fields:
            if "<{}>".format(field[0]) in line:
                temp = list(filter(None,line.split(" ")))
                retDict[field[0]] = field[1](temp[1])
    return(retDict)

def write_descriptor(   f_name,
                        lines,
                        J0):

    with open(f_name,'w') as f:
        for line in lines:
            if 'const T rest_fraction_D2Q5Constants<T>::J0 =' in line:
                f.write("{} {};\n".format(line,J0))
            else:
                f.write("{}\n".format(line))
lines = [ 
        # {{{
        '#include "rest_fraction_lattice.h"',
        'namespace plb {',
        'namespace descriptors {',
        '',
        '// AdvectionDiffusion D2Q5 //////////////////////////////////////////////',
        'template<typename T>',
        'const T rest_fraction_D2Q5Constants<T>::invD = (T)1 / (T) d;',
        '',
        'template<typename T>',
        'const int rest_fraction_D2Q5Constants<T>::vicinity = 1;',
        '',
        'template<typename T>',
        'const T rest_fraction_D2Q5Constants<T>::J0 =',
        '',
        'template<typename T>',
        'const int rest_fraction_D2Q5Constants<T>::c',
        '    [rest_fraction_D2Q5Constants<T>::q][rest_fraction_D2Q5Constants<T>::d] =',
        '{',
        '    { 0, 0},',
        '    {-1, 0},',
        '    {0, -1},',
        '    {1,0}, ',
        '    { 0,1}',
        '};',
        '',
        'template<typename T>',
        'const int rest_fraction_D2Q5Constants<T>::cNormSqr[rest_fraction_D2Q5Constants<T>::q] =',
        '{   0,',
        '    1,',
        '    1,',
        '    1,',
        '    1 ',
        '    };',
        '',
        'template<typename T>',
        'const T rest_fraction_D2Q5Constants<T>::t[rest_fraction_D2Q5Constants<T>::q] =',
        '    {',
        '        (T)J0,',
        '        (T)(1-J0)/4,',
        '        (T)(1-J0)/4,',
        '        (T)(1-J0)/4,',
        '        (T)(1-J0)/4',
        '    };',
        '',
        'template<typename T>',
        'const T rest_fraction_D2Q5Constants<T>::cs2 = (T)1 / (T)3;',
        '',
        'template<typename T>',
        'const T rest_fraction_D2Q5Constants<T>::invCs2 = (T)3;',
        '',
        'template<typename T>',
        'const char rest_fraction_Descriptor<T>::name[] = "rest_fraction_Descriptor";',
        '',
        '   }  // namespace descriptors',
        '}  // namespace plb'
        # }}}
        ]



fields = [
            ('K_0',float),
            ('J_0',float),
        ]
params = read_xml('params.xml',fields)

K_0 = params['K_0']

write_descriptor(   'rest_fraction_lattice.hh',
                    lines,
                    K_0)
