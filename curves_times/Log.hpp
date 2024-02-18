#include <iostream>
#include <string>
#include <memory>

class Log{
public:
    enum Level {Debug = 0, Info = 1, Warning = 2, Error = 3};

    static Level level;

    static void debug(const std::string& message){
        if(level == Debug)
            std::cout << "Debug: " << message << std::endl;
    }
    static void info(const std::string& message){
        if(level <= Info)
            std::cout << "Info: " << message << std::endl;
    }
    static void warning(const std::string& message){
        if(level <= Warning)
            std::cout << "Warning: " << message << std::endl;
    }
    static void error(const std::string& message){
        std::cout << "Error: " << message << std::endl;
        std::exit(1);
    }
};  Log::Level Log::level = Info;