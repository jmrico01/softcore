extern "C" int _fltused = 0;

extern "C" void* memset_rename_when_necessary(void* dest, int c, size_t count)
{
    unsigned char* d = (unsigned char*)dest;
    unsigned char value = *(char*)(&c);
    for (size_t i = 0; i < count; i++) {
        *(d++) = value;
    }
    
    return dest;
}