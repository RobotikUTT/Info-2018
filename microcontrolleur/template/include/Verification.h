#ifndef VERIFICATION_H
#define VERIFICATION_H

#include "../include/Interface.h" //This class need to know the class above



class Verification{
  private:
    Verification();


  public:
    static void CreateSingleton();
    ~Verification();
    void Update();
};
extern Verification* g_verification_instance_ptr;

#endif
