#include <iostream>
#include <string>
#include <vector>

using namespace std;

// typedef string mystringType;
typedef string mystringType;


class cartesian_product
{
 private:
  

  void cartesian_product_for_2strings(const vector<mystringType> &a, const vector<mystringType> &b, vector<mystringType> &res)
  {
    res.clear();
    for (unsigned i=0; i<a.size(); ++i)
    {
      for (unsigned j=0; j<b.size(); ++j)
      {
	mystringType tmp = a[i] + b[j];
	res.push_back(tmp);
      }
    }
  }

  void nth_cartesian_product_of_symbolset(unsigned replicates, const mystringType symbols, vector<mystringType> &res)
  {
    res.clear();

    vector<mystringType> start, tmp_res, tmp_res2;
    for (unsigned i=0; i < symbols.size(); ++i)
    {
      // string:
      mystringType tmp = mystringType(1u,symbols[i]);
      // faststring:
      //      mystringType tmp = mystringType(symbols[i], 1u);
      
      start.push_back(tmp);
    }
    
    //    cout << "DEBUG: start" << endl;
    //    print_container(cout, start, "", ",","");
    //    cout << "DEBUG: start -- END" << endl;
    
    tmp_res = start;
    for (unsigned j=1; j<replicates; ++j)
    {
      cartesian_product_for_2strings(tmp_res, start, tmp_res2);
      tmp_res = tmp_res2;
    }

    res = tmp_res2;
  }


 public:
  void operator()(unsigned replicates, const mystringType symbols, vector<mystringType> &res)
  {
    nth_cartesian_product_of_symbolset(replicates, symbols, res);
  }
};



