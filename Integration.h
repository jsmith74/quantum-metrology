#ifndef INTEGRATION_H_INCLUDED
#define INTEGRATION_H_INCLUDED

#include <vector>

class Integration{

    public:

        Integration();
        double integrate(std::vector<double>& f,int lowerBound,int higherBound);

    private:



};


#endif // INTEGRATION_H_INCLUDED
