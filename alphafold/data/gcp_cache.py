import hashlib
from google.cloud import storage
from alphafold.data import parsers
from alphafold.data.pipeline_multimer import _make_chain_id_map
import dataclasses
import json
import os
import shutil

def cache_hash(input: str) -> str:
    return hashlib.sha256(input.encode()).hexdigest()

def preload_multimer_cached_msas(input_fasta_path, msa_output_dir, cache_bucket='af2_cache', msa_gen_params="") -> None:    
    storage_client = storage.Client()
    bucket = storage_client.bucket(cache_bucket)

    with open(input_fasta_path) as f:
        input_fasta_str = f.read()
    input_seqs, input_descs = parsers.parse_fasta(input_fasta_str)

    chain_id_map = _make_chain_id_map(sequences=input_seqs,
                                        descriptions=input_descs)
    chain_id_map_path = os.path.join(msa_output_dir, 'chain_id_map.json')
    with open(chain_id_map_path, 'w') as f:
        chain_id_map_dict = {chain_id: dataclasses.asdict(fasta_chain)
                            for chain_id, fasta_chain in chain_id_map.items()}
        json.dump(chain_id_map_dict, f, indent=4, sort_keys=True)

    for chain_id, seq in chain_id_map.items():
        seq_hash = cache_hash(seq.sequence + msa_gen_params)
        remote_path = seq_hash + ".zip"
        blob = bucket.blob(remote_path)
        if blob.exists():
            archive_path = os.path.join(msa_output_dir, chain_id + ".zip")
            blob.download_to_filename(archive_path)
            local_path = os.path.join(msa_output_dir, chain_id)
            print(f"Loading precomputed MSA from {remote_path} as {local_path}")
            shutil.unpack_archive(archive_path, local_path)
            os.remove(archive_path)



def upload_multimer_zipped_msa(msa_output_dir, cache_bucket='af2_cache', msa_gen_params="") -> None:
    storage_client = storage.Client()
    bucket = storage_client.bucket(cache_bucket)

    chain_id_map_path = os.path.join(msa_output_dir, 'chain_id_map.json')
    with open(chain_id_map_path, 'r') as f:
        chain_id_map_dict = json.load(f)

    for chain_id, chain_data in chain_id_map_dict.items():
        seq_hash = cache_hash(chain_data['sequence'] + msa_gen_params) 
        remote_path = seq_hash + ".zip"
        blob = bucket.blob(remote_path)

        if not blob.exists():
            local_path = os.path.join(msa_output_dir, chain_id)
            archive_path = os.path.join(msa_output_dir, chain_id + ".zip")
            shutil.make_archive(local_path, 'zip', local_path)
            blob.upload_from_filename(archive_path)

            print(f"Uploaded zipped MSA from {local_path} to {remote_path}")
        else:
            print(f"File {remote_path} already exists in the cache, skipping upload.")

