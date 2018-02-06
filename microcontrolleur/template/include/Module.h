#ifndef MODULE_H
#define MODULE_H

#include "../include/Verification.h" //This class need to know the class above



class Module{
  private:
    Module();


  public:
    static void CreateSingleton();
    ~Module();
    void Update();
};
extern Module* g_module_instance_ptr;

#endif
