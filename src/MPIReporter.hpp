#pragma once

#include <string>
#include <vector>

#include "fwd.hpp"

// https://de.wikibooks.org/wiki/C%2B%2B-Programmierung:_Entwurfsmuster:_Singleton
class MPIReporter {
public:
    static MPIReporter* instance() {
        static CGuard g;  // Speicherbereinigung
        if (!_instance)
            _instance = new MPIReporter();
        return _instance;
    }
    void StoreMessage(int rank, std::string message) {
        messages.push_back("Message from processor: " + std::to_string(rank) + ": " + message);
    };
    std::vector<std::string> GetAllMessages() { return messages; };

private:
    std::vector<std::string> messages;
    static MPIReporter* _instance;
    MPIReporter() {}                 /* verhindert, dass ein Objekt von außerhalb von N erzeugt wird. */
                                     // protected, wenn man von der Klasse noch erben möchte
    MPIReporter(const MPIReporter&); /* verhindert, dass eine weitere Instanz via
Kopie-Konstruktor erstellt werden kann */
    ~MPIReporter() {}
    class CGuard {
    public:
        ~CGuard() {
            if (NULL != MPIReporter::_instance) {
                delete MPIReporter::_instance;
                MPIReporter::_instance = NULL;
            }
        }
    };
};