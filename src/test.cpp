// resizing vector
#include <iostream>
#include <vector>

int main ()
{
  std::vector<std::vector<int> > myvector(1, std::vector<int>(1));


  // set some initial content:
  for (int i=1;i<10;i++) myvector[0].push_back(i);

  myvector.resize(5);
  myvector.resize(4);
  //myvector.resize(8,100);
  //myvector.resize(12);

  std::cout << "myvector contains:";
  for (size_t i=0;i<myvector[0].size();i++)
    std::cout << ' ' << myvector[0][i];
  std::cout << '\n';

  return 0;
}
