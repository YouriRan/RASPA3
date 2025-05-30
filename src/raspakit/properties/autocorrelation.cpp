module;

#ifdef USE_LEGACY_HEADERS
#include <filesystem>
#include <fstream>
#include <iostream>
#include <source_location>
#include <sstream>
#include <string>
#include <vector>
#endif

module property_autocorrelation;

#ifndef USE_LEGACY_HEADERS
import <vector>;
import <string>;
import <filesystem>;
import <iostream>;
import <sstream>;
import <fstream>;
import <source_location>;
#endif

import archive;

PropertyAutocorrelation::PropertyAutocorrelation(size_t numberOfBuffersACF, size_t bufferLengthACF, size_t sampleEvery,
                                                 size_t writeEvery, std::string variableName)
    : numberOfBuffersACF(numberOfBuffersACF),
      bufferLengthACF(bufferLengthACF),
      sampleEvery(sampleEvery),
      writeEvery(writeEvery),
      variableName(variableName),
      buffer(numberOfBuffersACF, std::vector<double>(bufferLengthACF)),
      accumulatedAcf(bufferLengthACF),
      counts(numberOfBuffersACF)
{
  // trick to space the origins evenly (see for example Rapaport 2004)
  for (size_t currentBuffer = 0; currentBuffer < numberOfBuffersACF; ++currentBuffer)
  {
    counts[currentBuffer] = -static_cast<int64_t>(currentBuffer * bufferLengthACF / numberOfBuffersACF);
  }
}

void PropertyAutocorrelation::addSample(size_t currentCycle, double value)
{
  sumValue += value;
  if (currentCycle % sampleEvery != 0uz) return;

  for (size_t currentBuffer = 0; currentBuffer < numberOfBuffersACF; ++currentBuffer)
  {
    if (counts[currentBuffer] >= 0)
    {
      size_t index = static_cast<size_t>(counts[currentBuffer]);
      buffer[currentBuffer][index] = value;
    }
    counts[currentBuffer]++;
  }

  double mean = sumValue / static_cast<double>(currentCycle + 1);
  for (size_t currentBuffer = 0; currentBuffer < numberOfBuffersACF; ++currentBuffer)
  {
    if (counts[currentBuffer] == static_cast<std::make_signed_t<std::size_t>>(bufferLengthACF))
    {
      for (size_t t = 0; t < bufferLengthACF; ++t)
      {
        for (size_t i = 0; i < bufferLengthACF - t; ++i)
        {
          accumulatedAcf[t] += (buffer[currentBuffer][i] - mean) * (buffer[currentBuffer][i+t] - mean) / static_cast<double>(bufferLengthACF - t);
        }
      }
      counts[currentBuffer] = 0;
      ++countAccumulatedACF;
    }
  }
}

void PropertyAutocorrelation::writeOutput(size_t systemId, size_t currentCycle)
{
    if (currentCycle % writeEvery != 0uz) return;

    std::filesystem::create_directory("autocorrelation");

    std::ofstream stream_acf_output(std::format("autocorrelation/autocorrelation_{}.s{}.txt", variableName, systemId));

    stream_acf_output << std::format("# {}, number of counts: {}\n", variableName, countAccumulatedACF);
    stream_acf_output << "# column 1: cycle [-]\n";
    stream_acf_output << "# column 2: ck [-]\n";
    stream_acf_output << "# column 3: ck/c0 [-]\n";
    stream_acf_output << "# column 4: t_eff [-]\n";

    double fac = 1.0 / static_cast<double>(countAccumulatedACF);

    // start at -1 to compensate 2* c0/c0 without if statement
    double tEff = -1.0;
    for (size_t k = 0; k < bufferLengthACF; ++k)
    {
      if (k * sampleEvery < currentCycle)
      {
        double C = accumulatedAcf[k] / accumulatedAcf[0];
        tEff += 2.0 * C;
        stream_acf_output << std::format("{} {} {} {}\n", k * sampleEvery, fac * accumulatedAcf[k], C, static_cast<double>(sampleEvery) * tEff);
      }
    }
    stream_acf_output << std::endl;
}
