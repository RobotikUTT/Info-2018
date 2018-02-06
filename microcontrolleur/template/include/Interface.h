#ifndef INTERFACE_H
#define INTERFACE_H


class Interface{
  private:
    Interface();

  public:
    static void CreateSingleton();
    ~Interface();
    void Update();

};
extern Interface* interface_instance_ptr;

#endif
